import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass

from .case_matrix import (
    CaseContext,
    accumulate_species_counts_from_groups,
    build_case_label,
    discover_groups,
    fmt_num,
    normalize_force_field_sets,
    normalize_species_config,
    normalize_component_ratios,
    resolve_system_size,
)
from .config import load_config
from .mdp import modify_mdp
from .paths import get_output_root
from .slurm import (
    build_sbatch_directives,
    collect_slurm_status_map,
    get_mdrun_extra_args,
    get_slurm_step_launcher,
    parse_replica_job_ids,
    shell_join,
    submit_sbatch,
)
from .status_reporting import make_job_logger
from .system_setup import build_initial_box_if_needed, rebuild_replica_grompp_log

AVOGADRO = 6.02214076e23


@dataclass(frozen=True)
class ConcentrationOptimizerConfig:
    target_group: str
    target_molarity_mol_l: float
    tolerance_mol_l: float
    max_iterations: int
    reference_count: int
    initial_ratio_name: str
    force_field_name: str
    temperature: float
    output_subdir: str
    density_average_fraction: float
    max_weight_change_factor: float
    initial_density_guess_kg_m3: float | None
    density_analysis_tolerance: float
    density_analysis_show_points: int
    density_analysis_window_size: int
    density_analysis_mean_threshold: float
    density_analysis_distance: int
    density_analysis_min_consecutive_points: int
    density_analysis_time_crop: str | None


@dataclass(frozen=True)
class CompletedIteration:
    iteration_idx: int
    target_group_weight: float
    observed_molarity_mol_l: float
    observed_density_kg_m3: float | None
    error_mol_l: float
    within_tolerance: bool
    iter_root: str
    group_counts: dict
    component_ratio: dict
    selected_frame_time_ps: float | None
    selected_frame_volume_nm3: float


def load_optimizer_config(cfg):
    raw = cfg.get("concentration_optimizer")
    if not isinstance(raw, dict) or not raw.get("enabled", False):
        raise ValueError(
            "concentration_optimizer is not enabled in the config. "
            "Set [concentration_optimizer].enabled = true to use this tool."
        )
    optimizer_cfg = ConcentrationOptimizerConfig(
        target_group=str(raw["target_group"]),
        target_molarity_mol_l=float(raw["target_molarity_mol_l"]),
        tolerance_mol_l=float(raw["tolerance_mol_l"]),
        max_iterations=int(raw.get("max_iterations", 5)),
        reference_count=int(raw.get("reference_count", 150)),
        initial_ratio_name=str(raw["initial_ratio_name"]),
        force_field_name=str(raw["force_field_name"]),
        temperature=float(raw["temperature"]),
        output_subdir=str(raw.get("output_subdir", "concentration_optimizer")),
        density_average_fraction=float(raw.get("density_average_fraction", 0.2)),
        max_weight_change_factor=float(raw.get("max_weight_change_factor", 3.0)),
        initial_density_guess_kg_m3=(
            None if raw.get("initial_density_guess_kg_m3") in (None, "")
            else float(raw.get("initial_density_guess_kg_m3"))
        ),
        density_analysis_tolerance=float(raw.get("density_analysis_tolerance", 0.01)),
        density_analysis_show_points=int(raw.get("density_analysis_show_points", 20)),
        density_analysis_window_size=int(raw.get("density_analysis_window_size", 10)),
        density_analysis_mean_threshold=float(raw.get("density_analysis_mean_threshold", 0.1)),
        density_analysis_distance=int(raw.get("density_analysis_distance", 1)),
        density_analysis_min_consecutive_points=int(raw.get("density_analysis_min_consecutive_points", 1)),
        density_analysis_time_crop=(
            str(raw.get("density_analysis_time_crop")).strip()
            if raw.get("density_analysis_time_crop") not in (None, "")
            else None
        ),
    )
    if optimizer_cfg.target_molarity_mol_l <= 0:
        raise ValueError("concentration_optimizer.target_molarity_mol_l must be > 0")
    if optimizer_cfg.tolerance_mol_l <= 0:
        raise ValueError("concentration_optimizer.tolerance_mol_l must be > 0")
    if optimizer_cfg.max_iterations <= 0:
        raise ValueError("concentration_optimizer.max_iterations must be > 0")
    if optimizer_cfg.reference_count <= 0:
        raise ValueError("concentration_optimizer.reference_count must be > 0")
    if not 0.0 < optimizer_cfg.density_average_fraction <= 1.0:
        raise ValueError("concentration_optimizer.density_average_fraction must be in (0, 1]")
    if optimizer_cfg.max_weight_change_factor <= 1.0:
        raise ValueError("concentration_optimizer.max_weight_change_factor must be > 1")
    if optimizer_cfg.initial_density_guess_kg_m3 is not None and optimizer_cfg.initial_density_guess_kg_m3 <= 0:
        raise ValueError("concentration_optimizer.initial_density_guess_kg_m3 must be > 0 when provided")
    return optimizer_cfg


def resolve_optimizer_case(cfg, optimizer_cfg, group_keys, species_order):
    scr = cfg["screening"]
    ratio_entries = normalize_component_ratios(scr, group_keys)
    ff_entries = normalize_force_field_sets(scr, species_order)
    ratio_entry = next((entry for entry in ratio_entries if entry["name"] == optimizer_cfg.initial_ratio_name), None)
    if ratio_entry is None:
        raise KeyError(
            f"concentration_optimizer.initial_ratio_name='{optimizer_cfg.initial_ratio_name}' "
            "does not match any [[screening.component_ratios]] entry name"
        )
    ff_entry = next((entry for entry in ff_entries if entry["name"] == optimizer_cfg.force_field_name), None)
    if ff_entry is None:
        raise KeyError(
            f"concentration_optimizer.force_field_name='{optimizer_cfg.force_field_name}' "
            "does not match any [[screening.force_field]] entry name"
        )
    active_itps = ff_entry["species_itps"]
    label = build_case_label(
        optimizer_cfg.temperature,
        ratio_entry["name"],
        ff_entry["name"],
        ratio_entry["weights"],
        active_itps,
        species_order,
        group_keys,
    )
    return CaseContext(
        label=label,
        sys_path="",
        job_log_path="",
        ratio_entry=ratio_entry,
        temp=optimizer_cfg.temperature,
        ff_entry=ff_entry,
        active_itps=active_itps,
    )


def optimizer_root_dir(output_root, optimizer_cfg, case_ctx):
    opt_label = (
        f"concopt_T{fmt_num(optimizer_cfg.temperature)}"
        f"_{optimizer_cfg.initial_ratio_name}"
        f"_{optimizer_cfg.force_field_name}"
        f"_{optimizer_cfg.target_group}_{fmt_num(optimizer_cfg.target_molarity_mol_l)}M"
    )
    return os.path.join(output_root, optimizer_cfg.output_subdir, opt_label)


def iteration_dir(opt_root, iteration_idx):
    return os.path.join(opt_root, f"iter_{iteration_idx}")


def optimizer_metadata_path(iter_root):
    return os.path.join(iter_root, "iteration_metadata.json")


def completed_result_path(opt_root):
    return os.path.join(opt_root, "optimized_ratio.json")


def optimizer_summary_csv_path(opt_root):
    return os.path.join(opt_root, "concentration_optimization_history.csv")


def density_selection_metadata_path(iter_root):
    return os.path.join(iter_root, "3_npt", "density_selection.json")


def iteration_result_path(iter_root):
    return os.path.join(iter_root, "iteration_result.json")


def write_json(path, payload):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")


def read_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def parse_box_volume_nm3(gro_path):
    with open(gro_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    if not lines:
        raise ValueError(f"Could not read final GRO file '{gro_path}'")
    box_tokens = lines[-1].split()
    if len(box_tokens) == 3:
        a, b, c = [float(x) for x in box_tokens]
        return a * b * c
    if len(box_tokens) == 9:
        xx, yy, zz, xy, xz, yx, yz, zx, zy = [float(x) for x in box_tokens]
        ax, ay, az = xx, yx, zx
        bx, by, bz = xy, yy, zy
        cx, cy, cz = xz, yz, zz
        return abs(
            ax * (by * cz - bz * cy)
            - ay * (bx * cz - bz * cx)
            + az * (bx * cy - by * cx)
        )
    raise ValueError(f"Unsupported GRO box format in '{gro_path}': {lines[-1].strip()}")


def compute_group_molarity_mol_l(group_count, volume_nm3):
    if volume_nm3 <= 0:
        raise ValueError("Final volume must be > 0")
    volume_l = volume_nm3 * 1e-24
    return group_count / (AVOGADRO * volume_l)


def parse_itp_molar_mass(itp_path):
    section = None
    total_mass = 0.0
    found_atoms = False
    with open(itp_path, "r", encoding="utf-8", errors="ignore") as f:
        for raw_line in f:
            line = raw_line.split(";", 1)[0].strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                section = line.strip("[]").strip().lower()
                continue
            if section != "atoms":
                continue
            parts = line.split()
            if len(parts) < 8:
                raise ValueError(
                    f"ITP file '{itp_path}' has an [ atoms ] row without a readable mass column: '{raw_line.rstrip()}'"
                )
            try:
                total_mass += float(parts[7])
            except ValueError as exc:
                raise ValueError(
                    f"ITP file '{itp_path}' has a non-numeric mass value in [ atoms ]: '{raw_line.rstrip()}'"
                ) from exc
            found_atoms = True
    if not found_atoms or total_mass <= 0:
        raise ValueError(
            f"ITP file '{itp_path}' does not define a valid positive molecular mass in its [ atoms ] section"
        )
    return total_mass


def build_species_molar_masses(active_itps):
    masses = {}
    for species_key, itp_name in active_itps.items():
        masses[species_key] = parse_itp_molar_mass(os.path.join("inputs", itp_name))
    return masses


def build_group_molar_masses(group_defs, group_keys, species_molar_masses):
    group_masses = {}
    for group_key in group_keys:
        total = 0.0
        for species_key, coeff in group_defs[group_key]["stoichiometric_ratio"].items():
            total += float(coeff) * float(species_molar_masses[species_key])
        if total <= 0:
            raise ValueError(f"Computed non-positive molar mass for group '{group_key}'")
        group_masses[group_key] = total
    return group_masses


def read_density_from_edr(edr_path, average_fraction):
    if not os.path.exists(edr_path):
        return None
    try:
        import panedr
    except ModuleNotFoundError:
        return None
    try:
        df = panedr.edr_to_df(edr_path)
    except Exception:
        return None
    if "Density" not in df.columns or df.empty:
        return None
    density_series = df["Density"].dropna()
    if density_series.empty:
        return None
    tail_len = max(1, int(math.ceil(len(density_series) * average_fraction)))
    return float(density_series.tail(tail_len).mean())


def run_density_analysis(cfg, opt_cfg, iter_root):
    npt_dir = os.path.join(iter_root, "3_npt")
    edr_path = os.path.join(npt_dir, "npt.edr")
    if not os.path.exists(edr_path):
        return None
    analysis_script = os.path.join(
        os.getcwd(),
        cfg["project_settings"]["scripts_dir"],
        cfg["project_settings"]["density_analysis_script"],
    )
    metadata_path = density_selection_metadata_path(iter_root)
    cmd = [
        sys.executable,
        analysis_script,
        edr_path,
        "--num_replicas",
        "1",
        "--tolerance",
        str(opt_cfg.density_analysis_tolerance),
        "--show_points",
        str(opt_cfg.density_analysis_show_points),
        "--window_size",
        str(opt_cfg.density_analysis_window_size),
        "--mean_threshold",
        str(opt_cfg.density_analysis_mean_threshold),
        "--distance",
        str(opt_cfg.density_analysis_distance),
        "--min_consecutive_points",
        str(opt_cfg.density_analysis_min_consecutive_points),
        "--metadata-output",
        metadata_path,
    ]
    if opt_cfg.density_analysis_time_crop:
        cmd.extend(["--time_crop", opt_cfg.density_analysis_time_crop])
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode != 0 or not os.path.exists(metadata_path):
        print(
            f"WARNING: density analysis failed for '{iter_root}'. "
            "Falling back to the final NPT frame volume."
        )
        if res.stdout.strip():
            print(res.stdout.strip())
        if res.stderr.strip():
            print(res.stderr.strip())
        return None
    try:
        return read_json(metadata_path)
    except (OSError, json.JSONDecodeError):
        print(
            f"WARNING: density analysis metadata could not be parsed for '{iter_root}'. "
            "Falling back to the final NPT frame volume."
        )
        return None


def extract_selected_volume_from_density_analysis(cfg, opt_cfg, iter_root):
    metadata = run_density_analysis(cfg, opt_cfg, iter_root)
    if not metadata:
        return None
    selected_points = metadata.get("selected_points") or []
    if not selected_points:
        return None
    selected = selected_points[0]
    output_file = selected.get("output_file")
    if not output_file or not os.path.exists(output_file):
        return None
    return {
        "average_density": float(metadata["average_density"]),
        "selected_time": float(selected["time"]),
        "selected_volume_nm3": parse_box_volume_nm3(output_file),
    }


def estimate_initial_target_group_weight(opt_cfg, base_component_ratio, base_group_counts):
    current_weight = float(base_component_ratio[opt_cfg.target_group])
    if opt_cfg.initial_density_guess_kg_m3 is None:
        return current_weight
    scale_factor = float(opt_cfg.initial_density_guess_kg_m3) / 1000.0
    return max(1e-6, current_weight * scale_factor)


def build_density_driven_initial_guess(cfg, opt_cfg, group_defs, group_keys, order, active_itps, base_component_ratio):
    if opt_cfg.initial_density_guess_kg_m3 is None:
        return None

    target_group_count = int(opt_cfg.reference_count)
    target_group_moles = target_group_count / AVOGADRO
    volume_l = target_group_moles / float(opt_cfg.target_molarity_mol_l)
    volume_nm3 = volume_l * 1e24
    scale_factor = float(cfg["project_settings"].get("box_scale_factor", 1.0))
    box_size_nm = (scale_factor * volume_nm3) ** (1.0 / 3.0)
    density_g_l = float(opt_cfg.initial_density_guess_kg_m3)
    total_mass_g = density_g_l * volume_l
    if total_mass_g <= 0:
        raise ValueError("Initial density guess leads to non-positive total mass")

    species_molar_masses = build_species_molar_masses(active_itps)
    group_molar_masses = build_group_molar_masses(group_defs, group_keys, species_molar_masses)

    target_group_mass_g = (target_group_count / AVOGADRO) * group_molar_masses[opt_cfg.target_group]

    if target_group_mass_g >= total_mass_g:
        raise ValueError(
            "The requested target molarity, reference_count, and initial density guess are inconsistent: "
            "the target group alone exceeds the total mass implied by the guessed density."
        )

    other_group_keys = [key for key in group_keys if key != opt_cfg.target_group]
    other_weights = {key: float(base_component_ratio[key]) for key in other_group_keys}
    total_other_weight = sum(other_weights.values())
    remaining_mass_g = total_mass_g - target_group_mass_g

    group_counts = {key: 0 for key in group_keys}
    group_counts[opt_cfg.target_group] = target_group_count

    if other_group_keys:
        if total_other_weight <= 0:
            raise ValueError(
                "The selected initial ratio assigns zero total weight to all non-target groups, "
                "so the remaining mass implied by the guessed density cannot be distributed."
            )
        weighted_other_molar_mass = sum(
            other_weights[key] * group_molar_masses[key]
            for key in other_group_keys
        )
        if weighted_other_molar_mass <= 0:
            raise ValueError("Non-target groups produced a non-positive weighted molar mass")
        scale_moles = remaining_mass_g / weighted_other_molar_mass
        for key in other_group_keys:
            group_counts[key] = max(0, int(round(scale_moles * other_weights[key] * AVOGADRO)))

    species_counts = accumulate_species_counts_from_groups(group_counts, group_defs, order)
    component_ratio = {key: float(max(group_counts[key], 1 if key == opt_cfg.target_group else 0)) for key in group_keys}
    return {
        "box_size_nm": box_size_nm,
        "group_counts": group_counts,
        "species_counts": species_counts,
        "component_ratio": component_ratio,
        "species_molar_masses": species_molar_masses,
        "group_molar_masses": group_molar_masses,
        "initial_volume_l": volume_l,
        "initial_total_mass_g": total_mass_g,
        "reference_count": target_group_count,
    }


def write_optimizer_history(opt_root, completed_results):
    headers = [
        "iteration",
        "target_group_weight",
        "observed_molarity_mol_l",
        "error_mol_l",
        "within_tolerance",
        "observed_density_kg_m3",
        "selected_frame_time_ps",
        "selected_frame_volume_nm3",
        "iter_root",
    ]
    with open(optimizer_summary_csv_path(opt_root), "w", encoding="utf-8") as f:
        f.write(",".join(headers) + "\n")
        for result in completed_results:
            density_value = "" if result.observed_density_kg_m3 is None else str(result.observed_density_kg_m3)
            row = [
                str(result.iteration_idx),
                str(result.target_group_weight),
                str(result.observed_molarity_mol_l),
                str(result.error_mol_l),
                "yes" if result.within_tolerance else "no",
                density_value,
                "" if result.selected_frame_time_ps is None else str(result.selected_frame_time_ps),
                str(result.selected_frame_volume_nm3),
                result.iter_root,
            ]
            f.write(",".join(row) + "\n")


def write_iteration_result(result):
    payload = {
        "iteration_idx": result.iteration_idx,
        "target_group_weight": result.target_group_weight,
        "observed_molarity_mol_l": result.observed_molarity_mol_l,
        "observed_density_kg_m3": result.observed_density_kg_m3,
        "error_mol_l": result.error_mol_l,
        "within_tolerance": result.within_tolerance,
        "group_counts": result.group_counts,
        "component_ratio": result.component_ratio,
        "selected_frame_time_ps": result.selected_frame_time_ps,
        "selected_frame_volume_nm3": result.selected_frame_volume_nm3,
    }
    write_json(iteration_result_path(result.iter_root), payload)


def iteration_stage_inputs(iter_root, cfg, species_order, active_itps, species_cfg, counts, label, temp):
    topology_include = cfg["project_settings"].get("topology_forcefield_include", "forcefield.itp")
    stage_names = ["1_min", "2_press", "3_npt"]
    sim_cfg = cfg["simulation_settings"]
    press_cfg = sim_cfg.get("conc_opt_press", sim_cfg["press"])
    npt_cfg = sim_cfg["conc_opt_npt"]
    template_dir = cfg["project_settings"]["template_dir"]

    for stage in stage_names:
        dest = os.path.join(iter_root, stage)
        os.makedirs(dest, exist_ok=True)
        for filename in cfg["project_settings"]["global_ff_files"]:
            shutil.copy(f"inputs/{filename}", f"{dest}/{filename}")
        for key in species_order:
            shutil.copy(f"inputs/{active_itps[key]}", f"{dest}/{active_itps[key]}")
            shutil.copy(f"inputs/{species_cfg[key]['gro']}", f"{dest}/{species_cfg[key]['gro']}")

        with open(os.path.join(dest, "topol.top"), "w", encoding="utf-8") as top:
            top.write(f'#include "{topology_include}"\n')
            for key in species_order:
                top.write(f'#include "{active_itps[key]}"\n')
            top.write(f"\n[ system ]\n{label}\n\n[ molecules ]\n")
            for key in species_order:
                top.write(f"{species_cfg[key]['resname']:<15} {counts[key]}\n")

        if stage == "1_min":
            shutil.copy(f"{template_dir}/min.mdp", f"{dest}/min.mdp")
        elif stage == "2_press":
            modify_mdp(
                f"{template_dir}/press.mdp",
                f"{dest}/press.mdp",
                {
                    "ref_t": press_cfg["ref_t"],
                    "nsteps": press_cfg["nsteps"],
                    "dt": press_cfg["dt"],
                    "ref_p": press_cfg["ref_p"],
                },
            )
        elif stage == "3_npt":
            modify_mdp(
                f"{template_dir}/npt.mdp",
                f"{dest}/npt.mdp",
                {
                    "ref_t": npt_cfg.get("ref_t", temp),
                    "nsteps": npt_cfg["nsteps"],
                    "dt": npt_cfg["dt"],
                    "ref_p": npt_cfg["ref_p"],
                },
            )


def write_optimizer_iteration_sh(iter_root, cfg, aggregate_log_path, config_path):
    project_cfg, slurm_cfg = cfg["project_settings"], cfg["slurm_settings"]
    module_block = "\n".join(project_cfg.get("gmx_modules", []))
    optimizer_entrypoint = os.path.join(os.getcwd(), "optimize_concentration.py")
    grompp_launcher = get_slurm_step_launcher(slurm_cfg, "setup", "grompp")
    mdrun_launcher = get_slurm_step_launcher(slurm_cfg, "setup", "mdrun")
    mdrun_extra_args = get_mdrun_extra_args(slurm_cfg, "setup")
    sbatch_block = build_sbatch_directives(
        slurm_cfg,
        "setup",
        f"COPT_{project_cfg['system_name']}",
        f"{iter_root}/concentration.out",
        f"{iter_root}/concentration.err",
    )
    script_content = f"""#!/bin/bash
{sbatch_block}

set -euo pipefail

{module_block}
source ~/miniconda3/bin/activate {project_cfg['conda_env']}

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

BASE="{iter_root}"
GMX="{project_cfg.get('gmx_executable_path', 'gmx_mpi')}"
AGG_LOG="{aggregate_log_path}"

append_stage_log() {{
    local stage_label="$1"
    local file_path="$2"
    if [ -f "$file_path" ]; then
        {{
            printf "\\n====================================================================================================\\n"
            printf "%s\\n" "$stage_label"
            printf "Source: %s\\n" "$file_path"
            printf "====================================================================================================\\n"
            cat "$file_path"
            printf "\\n"
        }} >> "$AGG_LOG"
    fi
}}

cd $BASE/1_min
if ! {shell_join([grompp_launcher, "$GMX", "grompp -f min.mdp -c start.gro -p topol.top -o min.tpr"])} > $BASE/1_min/grompp_min.log 2>&1; then
    append_stage_log "$(basename "$BASE") | 1_min/grompp_min.log" "$BASE/1_min/grompp_min.log"
    echo "ERROR: concentration optimizer grompp failed for minimization."
    exit 1
fi
append_stage_log "$(basename "$BASE") | 1_min/grompp_min.log" "$BASE/1_min/grompp_min.log"
{shell_join([mdrun_launcher, "$GMX", "mdrun -deffnm min -c min_out.gro", mdrun_extra_args])}

cd $BASE/2_press
if ! {shell_join([grompp_launcher, "$GMX", "grompp -f press.mdp -c ../1_min/min_out.gro -p topol.top -o press.tpr"])} > $BASE/2_press/grompp_press.log 2>&1; then
    append_stage_log "$(basename "$BASE") | 2_press/grompp_press.log" "$BASE/2_press/grompp_press.log"
    echo "ERROR: concentration optimizer grompp failed for pressure equilibration."
    exit 1
fi
append_stage_log "$(basename "$BASE") | 2_press/grompp_press.log" "$BASE/2_press/grompp_press.log"
{shell_join([mdrun_launcher, "$GMX", "mdrun -deffnm press -c press_out.gro -rdd 1.2 -pin on", mdrun_extra_args])}

cd $BASE/3_npt
if ! {shell_join([grompp_launcher, "$GMX", "grompp -f npt.mdp -c ../2_press/press_out.gro -p topol.top -o npt.tpr"])} > $BASE/3_npt/grompp_npt.log 2>&1; then
    append_stage_log "$(basename "$BASE") | 3_npt/grompp_npt.log" "$BASE/3_npt/grompp_npt.log"
    echo "ERROR: concentration optimizer grompp failed for NPT."
    exit 1
fi
append_stage_log "$(basename "$BASE") | 3_npt/grompp_npt.log" "$BASE/3_npt/grompp_npt.log"
{shell_join([mdrun_launcher, "$GMX", "mdrun -deffnm npt -c npt_out.gro -pin on", mdrun_extra_args])}

python3 "{optimizer_entrypoint}" "{config_path}" --auto-continue --dependency-job-id "$SLURM_JOB_ID"
"""
    with open(os.path.join(iter_root, "run_concentration_iteration.sh"), "w", encoding="utf-8") as f:
        f.write(script_content)


def completed_iteration_indices(opt_root):
    indices = []
    if not os.path.exists(opt_root):
        return indices
    for name in os.listdir(opt_root):
        if name.startswith("iter_"):
            suffix = name.split("_", 1)[1]
            if suffix.isdigit():
                indices.append(int(suffix))
    return sorted(indices)


def iteration_job_status(iter_root):
    log_path = os.path.join(iter_root, "launch.log")
    job_ids = parse_replica_job_ids(log_path, 1)
    job_id = job_ids["setup_job_id"]
    if not job_id:
        return None, "not_submitted"
    status = collect_slurm_status_map([job_id]).get(job_id, "unknown")
    return job_id, status


def load_completed_iteration(cfg, opt_cfg, iter_root):
    metadata = read_json(optimizer_metadata_path(iter_root))
    npt_out_gro = os.path.join(iter_root, "3_npt", "npt_out.gro")
    if not os.path.exists(npt_out_gro):
        return None
    density_selection = extract_selected_volume_from_density_analysis(cfg, opt_cfg, iter_root)
    if density_selection:
        volume_nm3 = density_selection["selected_volume_nm3"]
        density = density_selection["average_density"]
        selected_time = density_selection["selected_time"]
    else:
        volume_nm3 = parse_box_volume_nm3(npt_out_gro)
        density = read_density_from_edr(
            os.path.join(iter_root, "3_npt", "npt.edr"),
            opt_cfg.density_average_fraction,
        )
        selected_time = None
    observed_molarity = compute_group_molarity_mol_l(metadata["group_counts"][opt_cfg.target_group], volume_nm3)
    error = observed_molarity - opt_cfg.target_molarity_mol_l
    return CompletedIteration(
        iteration_idx=int(metadata["iteration_idx"]),
        target_group_weight=float(metadata["target_group_weight"]),
        observed_molarity_mol_l=observed_molarity,
        observed_density_kg_m3=density,
        error_mol_l=error,
        within_tolerance=abs(error) <= opt_cfg.tolerance_mol_l,
        iter_root=iter_root,
        group_counts=metadata["group_counts"],
        component_ratio=metadata["component_ratio"],
        selected_frame_time_ps=selected_time,
        selected_frame_volume_nm3=volume_nm3,
    )


def next_weight_guess(opt_cfg, completed_results):
    if not completed_results:
        raise ValueError("At least one completed result is required to compute the next weight guess")
    latest = completed_results[-1]
    target = opt_cfg.target_molarity_mol_l
    if latest.observed_molarity_mol_l <= 0:
        raise ValueError("Observed molarity must be > 0 to compute the next weight guess")
    if len(completed_results) == 1:
        raw_guess = latest.target_group_weight * (target / latest.observed_molarity_mol_l)
    else:
        prev = completed_results[-2]
        denom = latest.observed_molarity_mol_l - prev.observed_molarity_mol_l
        if abs(denom) < 1e-12:
            raw_guess = latest.target_group_weight * (target / latest.observed_molarity_mol_l)
        else:
            raw_guess = latest.target_group_weight + (
                (target - latest.observed_molarity_mol_l)
                * (latest.target_group_weight - prev.target_group_weight)
                / denom
            )
    max_factor = opt_cfg.max_weight_change_factor
    lower = latest.target_group_weight / max_factor
    upper = latest.target_group_weight * max_factor
    bounded = min(max(raw_guess, lower), upper)
    return max(1e-6, bounded)


def render_ratio_snippet(component_ratio):
    lines = ["[[screening.component_ratios]]"]
    for key, value in component_ratio.items():
        lines.append(f'{key} = {value:.8g}')
    return "\n".join(lines)


def prepare_iteration(
    cfg,
    opt_cfg,
    opt_root,
    iter_idx,
    case_ctx,
    species_cfg,
    order,
    group_defs,
    group_keys,
    component_ratio,
    gmx_path,
    config_path,
    dependency_job_id=None,
    explicit_species_counts=None,
    explicit_group_counts=None,
    explicit_box_size_nm=None,
):
    iter_root = iteration_dir(opt_root, iter_idx)
    os.makedirs(iter_root, exist_ok=True)
    job_log, log_path = make_job_logger(iter_root, f"CONCENTRATION OPTIMIZER ITERATION {iter_idx}")
    grompp_log_path = rebuild_replica_grompp_log(iter_root)

    if explicit_species_counts is not None and explicit_group_counts is not None and explicit_box_size_nm is not None:
        counts = explicit_species_counts
        sizing_info = {
            "mode": "concentration_optimizer_density_guess",
            "n_total_components": sum(int(v) for v in explicit_group_counts.values()),
            "total_atoms": None,
            "box_size_nm": float(explicit_box_size_nm),
            "box_source": "concentration_optimizer.initial_density_guess_kg_m3",
            "component_counts": explicit_group_counts,
        }
    else:
        counts, sizing_info = resolve_system_size(cfg, species_cfg, order, group_defs, group_keys, component_ratio)
    label = (
        f"concopt_iter_{iter_idx}_"
        f"T{fmt_num(case_ctx.temp)}_{case_ctx.ff_entry['name']}_{opt_cfg.target_group}_{fmt_num(component_ratio[opt_cfg.target_group])}"
    )

    iteration_stage_inputs(iter_root, cfg, order, case_ctx.active_itps, species_cfg, counts, label, case_ctx.temp)
    build_initial_box_if_needed(
        iter_root,
        label,
        counts,
        order,
        species_cfg,
        gmx_path,
        sizing_info["box_size_nm"],
        job_log,
        grompp_log_path,
        iter_idx,
    )

    metadata = {
        "iteration_idx": iter_idx,
        "target_group": opt_cfg.target_group,
        "target_group_weight": component_ratio[opt_cfg.target_group],
        "component_ratio": component_ratio,
        "group_counts": sizing_info["component_counts"],
        "species_counts": counts,
        "box_size_nm": sizing_info["box_size_nm"],
        "force_field_name": case_ctx.ff_entry["name"],
        "temperature": case_ctx.temp,
    }
    write_json(optimizer_metadata_path(iter_root), metadata)
    write_optimizer_iteration_sh(iter_root, cfg, grompp_log_path, config_path)
    try:
        job_id = submit_sbatch(
            os.path.join(iter_root, "run_concentration_iteration.sh"),
            dependency_job_id=dependency_job_id,
        )
    except RuntimeError as exc:
        job_log(f"Iteration {iter_idx}: ERROR submitting concentration optimization job: {exc}")
        print(f"ERROR submitting concentration optimization iteration {iter_idx}: {exc}")
        return 1
    job_log(f"Replica 1: submitted setup job with id {job_id}")
    print(
        f"Submitted concentration optimization iteration {iter_idx} with job id {job_id} "
        f"(target-group weight={component_ratio[opt_cfg.target_group]:.8g})"
    )
    return 0


def collect_completed_results(cfg, opt_cfg, opt_root):
    results = []
    for iter_idx in completed_iteration_indices(opt_root):
        iter_root = iteration_dir(opt_root, iter_idx)
        result = load_completed_iteration(cfg, opt_cfg, iter_root)
        if result is not None:
            results.append(result)
    return results


def print_completed_summary(opt_cfg, completed_results):
    if not completed_results:
        print("No completed concentration-optimization iterations found yet.")
        return
    print("Completed concentration-optimization iterations:")
    for result in completed_results:
        density_txt = "n/a" if result.observed_density_kg_m3 is None else f"{result.observed_density_kg_m3:.3f} kg/m^3"
        print(
            f"- iter_{result.iteration_idx}: "
            f"weight={result.target_group_weight:.8g}, "
                f"molarity={result.observed_molarity_mol_l:.6f} mol/L, "
                f"error={result.error_mol_l:+.6f} mol/L, "
                f"density={density_txt}, "
                f"selected_frame_time_ps={'n/a' if result.selected_frame_time_ps is None else f'{result.selected_frame_time_ps:.3f}'}, "
                f"within_tolerance={'yes' if result.within_tolerance else 'no'}"
            )


def run_optimizer(config_path, auto_continue=False, dependency_job_id=None):
    cfg = load_config(config_path)
    opt_cfg = load_optimizer_config(cfg)
    species_cfg = normalize_species_config(cfg["species"])
    group_defs, group_keys, active_species_keys = discover_groups(cfg, species_cfg)
    if opt_cfg.target_group not in group_keys:
        raise KeyError(
            f"concentration_optimizer.target_group='{opt_cfg.target_group}' is not defined under [groups.*]"
        )

    order = active_species_keys
    case_ctx = resolve_optimizer_case(cfg, opt_cfg, group_keys, order)
    output_root = get_output_root(cfg)
    opt_root = optimizer_root_dir(output_root, opt_cfg, case_ctx)
    os.makedirs(opt_root, exist_ok=True)

    completed_results = collect_completed_results(cfg, opt_cfg, opt_root)
    write_optimizer_history(opt_root, completed_results)
    print_completed_summary(opt_cfg, completed_results)
    for result in completed_results:
        write_iteration_result(result)

    if completed_results and completed_results[-1].within_tolerance:
        converged = completed_results[-1]
        payload = {
            "target_group": opt_cfg.target_group,
            "target_molarity_mol_l": opt_cfg.target_molarity_mol_l,
            "tolerance_mol_l": opt_cfg.tolerance_mol_l,
            "converged_iteration": converged.iteration_idx,
            "observed_molarity_mol_l": converged.observed_molarity_mol_l,
            "observed_density_kg_m3": converged.observed_density_kg_m3,
            "selected_frame_time_ps": converged.selected_frame_time_ps,
            "selected_frame_volume_nm3": converged.selected_frame_volume_nm3,
            "component_ratio": converged.component_ratio,
            "group_counts": converged.group_counts,
            "screening_component_ratio_snippet": render_ratio_snippet(converged.component_ratio),
        }
        write_json(completed_result_path(opt_root), payload)
        print("\nTarget molarity reached within tolerance.")
        print(f"Optimized ratio file written to {completed_result_path(opt_root)}")
        print("\nSuggested [[screening.component_ratios]] snippet:\n")
        print(render_ratio_snippet(converged.component_ratio))
        return 0

    if len(completed_results) >= opt_cfg.max_iterations:
        print(
            f"\nMaximum number of iterations reached ({opt_cfg.max_iterations}) without hitting the target tolerance."
        )
        if completed_results:
            best = min(completed_results, key=lambda item: abs(item.error_mol_l))
            print(
                f"Closest result was iter_{best.iteration_idx}: "
                f"{best.observed_molarity_mol_l:.6f} mol/L "
                f"(error {best.error_mol_l:+.6f} mol/L)."
            )
            print("\nClosest [[screening.component_ratios]] snippet:\n")
            print(render_ratio_snippet(best.component_ratio))
        return 0 if auto_continue else 1

    existing_indices = completed_iteration_indices(opt_root)
    if existing_indices:
        latest_idx = max(existing_indices)
        latest_iter_root = iteration_dir(opt_root, latest_idx)
        latest_result = load_completed_iteration(cfg, opt_cfg, latest_iter_root)
        if latest_result is None:
            job_id, status = iteration_job_status(latest_iter_root)
            print(
                f"\nLatest iteration iter_{latest_idx} has not completed yet. "
                f"Slurm job id: {job_id or 'n/a'}, status: {status}."
            )
            return 0

    gmx_path = cfg["project_settings"].get("gmx_executable_path")
    if not gmx_path or not os.path.exists(gmx_path):
        gmx_path = shutil.which("gmx_mpi")
    if not gmx_path:
        print("ERROR: GMX not found")
        return 1

    if not completed_results:
        next_iter_idx = 1
        base_component_ratio = dict(case_ctx.ratio_entry["weights"])
        density_driven_guess = build_density_driven_initial_guess(
            cfg,
            opt_cfg,
            group_defs,
            group_keys,
            order,
            case_ctx.active_itps,
            base_component_ratio,
        )
        if density_driven_guess is not None:
            component_ratio = density_driven_guess["component_ratio"]
        else:
            component_ratio = dict(base_component_ratio)
            counts, sizing_info = resolve_system_size(cfg, species_cfg, order, group_defs, group_keys, component_ratio)
            component_ratio[opt_cfg.target_group] = estimate_initial_target_group_weight(
                opt_cfg,
                component_ratio,
                sizing_info["component_counts"],
            )
    else:
        next_iter_idx = completed_results[-1].iteration_idx + 1
        component_ratio = dict(completed_results[-1].component_ratio)
        component_ratio[opt_cfg.target_group] = next_weight_guess(opt_cfg, completed_results)

    print(
        f"\nPreparing concentration-optimization iteration {next_iter_idx} "
        f"for target group '{opt_cfg.target_group}' at {opt_cfg.target_molarity_mol_l:.6f} mol/L."
    )
    return prepare_iteration(
        cfg,
        opt_cfg,
        opt_root,
        next_iter_idx,
        case_ctx,
        species_cfg,
        order,
        group_defs,
        group_keys,
        component_ratio,
        gmx_path,
        os.path.abspath(config_path),
        dependency_job_id=dependency_job_id,
        explicit_species_counts=(density_driven_guess["species_counts"] if not completed_results and density_driven_guess is not None else None),
        explicit_group_counts=(density_driven_guess["group_counts"] if not completed_results and density_driven_guess is not None else None),
        explicit_box_size_nm=(density_driven_guess["box_size_nm"] if not completed_results and density_driven_guess is not None else None),
    )


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", nargs="?", default="config.toml")
    parser.add_argument("--auto-continue", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--dependency-job-id", help=argparse.SUPPRESS)
    args = parser.parse_args()
    sys.exit(
        run_optimizer(
            args.config,
            auto_continue=args.auto_continue,
            dependency_job_id=args.dependency_job_id,
        )
    )
