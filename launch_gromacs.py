import os, shutil, subprocess, itertools, sys, argparse, re, math, csv
from datetime import datetime

try:
    import tomllib
except ModuleNotFoundError:
    try:
        import tomli as tomllib
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "Python < 3.11 requires 'tomli' for TOML config parsing. Install it with: pip install tomli"
        ) from exc

#This script 
#
#


def create_default_config(filename):
    """Creates a template configuration TOML file if it doesn't exist."""
    default_cfg = """# GROMACS launch configuration (TOML)
# Define every molecular species using [[species]] entries.
# The script derives `resname` and `atoms_per_molecule` automatically from each species GRO file.
# Define formulation groups under [groups.*] using species names and stoichiometric coefficients.
# System molar ratios are specified only between groups.
#
# Screened system compositions are defined in [[screening.component_ratios]]
# as weights over groups. Group internal stoichiometry stays in [groups.*].stoichiometric_ratio.
# Force-field choices (species -> ITP file mappings) are screened independently in
# [[screening.force_field]], allowing a Cartesian product of ratios x temperatures x FF sets.

[project_settings]
system_name = "LiFSI_DME_TOL"
num_replicas = 1
scripts_dir = "scripts"
density_analysis_script = "density_analysis_tool.py"
density_analysis_time_crop = "600:1000"
box_size_nm = 15.0
box_scale_factor = 1.5
num_prod_chunks = 4
gmx_modules = [
#    "module purge",
#    "module load bsc/1.0",
#    "module load intel/2024.2",
#    "module load impi/2021.13",
#    "module load mkl/2024.2",
    "module load oneapi/2024.2",
    "module load gromacs/2024.2"
]
conda_env = "mdanalysis"
global_ff_files = ["ffbonded.itp", "ffnonbonded.itp", "forcefield.itp"]
template_dir = "mdp_templates"
gmx_executable_path = "/apps/GPP/GROMACS/2024.2/INTEL24.2/IMPI21.13/bin/gmx_mpi"

[system_sizing]
# Options:
#   "target_atoms"      -> estimate counts from target atom count + ratios
#   "reference_component" -> keep a chosen top-level component count fixed (salt or solvent)

#mode = "target_atoms"
#target_atoms = 12000
#estimate_box_from_atoms = true # -> if set to true, insert-molecules will use this volume
#atom_number_density_atoms_per_nm3 = 85.0

mode = "reference_component"
# Keep one top-level component count fixed using the same names as in
# [groups.*] (example names are arbitrary):
reference_component_key = "LiFSI_salt"    # or "ZnTFSI2_salt", "Solvents", etc.
reference_component_count = 150
estimate_box_from_atoms = false # -> if set to true, insert-molecules will use this volume
atom_number_density_atoms_per_nm3 = 85.0
#
# If estimate_box_from_atoms = false, project_settings.box_size_nm is used.

[screening]
# Component ratio basis = arbitrary weights across groups.
# Example below corresponds to LiFSI_salt:DME_solv:TOL_solv = 1:2:2.
target_temps = [298.15]

[[screening.component_ratios]]
# Weights over groups, not necessarily normalized.
LiFSI_salt = 1.0
DME_solv = 2.0
TOL_solv = 2.0

# Add more entries to screen more total formulations (example with grouped solvents):
# [[screening.component_ratios]]
# LiFSI_salt = 1.0
# Solvents = 1.0

[[screening.force_field]]
name = "ff_set1"
cation = "li.itp"
anion = "fsi.itp"
solvent_1 = "dme_qforce_cm5_1.itp"
solvent_2 = "tol_qforce_cm5_1.itp"

[[screening.force_field]]
name = "ff_set2"
cation = "li.itp"
anion = "fsi.itp"
solvent_1 = "dme_qforce_cm5_09.itp"
solvent_2 = "tol_qforce_cm5_09.itp"

[slurm_settings]
account = "cice86"
partition_setup = "gpp"
partition_prod = "gpp"
qos_setup = "gp_debug"
qos_prod = "gp_resa"
time_setup = "02:00:00"
time_prod = "72:00:00"
nodes_setup = 1
ntasks_per_node_setup = 56
cpus_per_task_setup = 2
nodes_prod = 1
ntasks_per_node_prod = 56
cpus_per_task_prod = 2

[[species]]
# Species names are arbitrary. Group stoichiometric_ratio references these names.
name = "cation"
gro = "li.gro"

[[species]]
name = "anion"
gro = "fsi_freq.gro"

[[species]]
name = "solvent_1"
gro = "dme.gro"

[[species]]
name = "solvent_2"
gro = "tol.gro"

# Add more species as needed (example):
# [[species]]
# name = "acn"
# gro = "acn.gro"

[groups.LiFSI_salt]
# One unit of the LiFSI salt group contains:
stoichiometric_ratio = { cation = 1, anion = 1 }

[groups.DME_solv]
stoichiometric_ratio = { solvent_1 = 1 }

[groups.TOL_solv]
stoichiometric_ratio = { solvent_2 = 1 }

# Alternative grouped solvent blend example:
# [groups.Solvents]
# stoichiometric_ratio = { solvent_1 = 2, solvent_2 = 1 }

[simulation_settings.press]
ref_p = 10.0
ref_t = 10.0
nsteps = 150000
dt = 0.002

[simulation_settings.anneal]
ref_p = 1.0
temp_high = 500.0
annealing_time = "0 200 400 1000"
nsteps = 500000
dt = 0.002

[simulation_settings.prod]
nsteps = 500000
dt = 0.002
"""
    with open(filename, "w") as f:
        f.write(default_cfg)
    print(f"\n[!] Template '{filename}' created.")
    sys.exit(0)

def modify_mdp(template_path, output_path, replacements):
    """Replaces MDP keys with tailored values."""
    with open(template_path, 'r') as f: lines = f.readlines()
    new_lines = []
    for line in lines:
        cleaned = line.split(';')[0].strip()
        if "=" in cleaned:
            key = cleaned.split('=')[0].strip()
            if key in replacements:
                new_lines.append(f"{key:<24}= {replacements[key]} ; Tailored\n")
                continue
        new_lines.append(line)
    with open(output_path, 'w') as f: f.writelines(new_lines)

def read_gro_species_metadata(gro_path):
    with open(gro_path, "r") as f:
        lines = f.readlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid GRO file '{gro_path}': too few lines")
    try:
        atom_count = int(lines[1].strip())
    except ValueError as exc:
        raise ValueError(f"Invalid GRO file '{gro_path}': line 2 must be atom count") from exc
    atom_lines = lines[2:2 + atom_count]
    if len(atom_lines) != atom_count:
        raise ValueError(f"Invalid GRO file '{gro_path}': expected {atom_count} atom lines")
    residue_names = {line[5:10].strip() for line in atom_lines if len(line) >= 10}
    if not residue_names:
        raise ValueError(f"Invalid GRO file '{gro_path}': could not parse residue name")
    if len(residue_names) != 1:
        raise ValueError(
            f"GRO file '{gro_path}' contains multiple residue names ({', '.join(sorted(residue_names))}); "
            "cannot auto-derive a single species name for topology output"
        )
    return {"resname": next(iter(residue_names)), "atoms_per_molecule": atom_count}

def normalize_species_config(raw_species, inputs_dir="inputs"):
    if not isinstance(raw_species, list) or not raw_species:
        raise ValueError(
            "species must be defined as a non-empty array of tables using [[species]] entries"
        )
    species_cfg = {}
    for entry in raw_species:
        if not isinstance(entry, dict):
            raise ValueError("Each [[species]] entry must be a table")
        name = entry.get("name")
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Each [[species]] entry must define a non-empty string 'name'")
        key = name.strip()
        if key in species_cfg:
            raise ValueError(f"Duplicate species name '{key}' in [[species]] entries")
        normalized = dict(entry)
        normalized["name"] = key
        gro_name = normalized.get("gro")
        if not isinstance(gro_name, str) or not gro_name.strip():
            raise ValueError(f"Species '{key}' must define a non-empty 'gro' filename")
        gro_name = gro_name.strip()
        normalized["gro"] = gro_name
        gro_meta = read_gro_species_metadata(os.path.join(inputs_dir, gro_name))
        normalized["resname"] = gro_meta["resname"]
        normalized["atoms_per_molecule"] = gro_meta["atoms_per_molecule"]
        species_cfg[key] = normalized
    return species_cfg

def discover_groups(cfg, species_cfg):
    groups_cfg = cfg.get("groups")
    if not groups_cfg or not isinstance(groups_cfg, dict):
        raise ValueError("[groups] must be defined as a table of named groups")
    group_defs = {}
    active_species_keys = []
    for group_key, group_entry in groups_cfg.items():
        if not isinstance(group_entry, dict):
            raise ValueError(f"groups.{group_key} must be a table")
        stoich = group_entry.get("stoichiometric_ratio")
        if not isinstance(stoich, dict) or not stoich:
            raise ValueError(f"groups.{group_key}.stoichiometric_ratio must be a non-empty table")
        norm = {}
        for sp_key, coeff in stoich.items():
            if sp_key not in species_cfg:
                raise KeyError(f"groups.{group_key}.stoichiometric_ratio references unknown species '{sp_key}'")
            coeff_i = int(coeff)
            if coeff_i <= 0:
                raise ValueError(f"groups.{group_key}.stoichiometric_ratio.{sp_key} must be > 0")
            norm[sp_key] = coeff_i
            if sp_key not in active_species_keys:
                active_species_keys.append(sp_key)
        group_defs[group_key] = {"stoichiometric_ratio": norm}
    if not active_species_keys:
        raise ValueError("[groups] must reference at least one species")
    return group_defs, list(group_defs.keys()), active_species_keys

def normalize_component_ratios(scr, component_keys):
    comps = scr.get("component_ratios")
    if not isinstance(comps, list) or not comps:
        raise KeyError(
            "Missing screening.component_ratios. "
            "Use [[screening.component_ratios]] with weights for groups."
        )
    out = []
    for idx, comp in enumerate(comps, start=1):
        if not isinstance(comp, dict):
            raise ValueError("Each screening.component_ratios entry must be a table")
        entry_name = comp.get("name", f"ratio_{idx}")
        unknown = [k for k in comp if k not in component_keys and k != "name"]
        if unknown:
            raise KeyError(f"Unknown component keys in ratio: {', '.join(unknown)}")
        weights = {k: float(comp.get(k, 0.0)) for k in component_keys}
        if sum(weights.values()) <= 0:
            raise ValueError("Each component ratio entry must have positive total weight")
        out.append({"name": str(entry_name), "weights": weights})
    return out

def normalize_force_field_sets(scr, species_keys):
    ff_sets = scr.get("force_field")
    if not isinstance(ff_sets, list) or not ff_sets:
        raise KeyError(
            "Missing screening.force_field. "
            "Use [[screening.force_field]] entries mapping species names to ITP files."
        )
    out = []
    for idx, entry in enumerate(ff_sets, start=1):
        if not isinstance(entry, dict):
            raise ValueError("Each screening.force_field entry must be a table")
        name = str(entry.get("name", f"ff_{idx}"))
        unknown = [k for k in entry if k not in species_keys and k != "name"]
        if unknown:
            raise KeyError(f"Unknown species keys in screening.force_field '{name}': {', '.join(unknown)}")
        missing = [k for k in species_keys if k not in entry]
        if missing:
            raise KeyError(f"Missing species ITP mappings in screening.force_field '{name}': {', '.join(missing)}")
        species_itps = {}
        for k in species_keys:
            v = entry[k]
            if not isinstance(v, str) or not v.strip():
                raise ValueError(f"screening.force_field '{name}' entry for species '{k}' must be a non-empty string")
            species_itps[k] = v.strip()
        out.append({"name": name, "species_itps": species_itps})
    return out

def allocate_weighted_counts(n_total, composition):
    total_weight = sum(float(v) for v in composition.values())
    if total_weight <= 0:
        raise ValueError("Composition total weight must be positive")
    raw = {k: (n_total * float(v) / total_weight) for k, v in composition.items()}
    counts = {k: int(v) for k, v in raw.items()}
    remainder = n_total - sum(counts.values())
    if remainder > 0:
        rank = sorted(raw, key=lambda k: (raw[k] - counts[k]), reverse=True)
        for k in rank[:remainder]:
            counts[k] += 1
    return counts

def accumulate_species_counts_from_groups(group_counts, group_defs, active_species_keys):
    species_counts = {k: 0 for k in active_species_keys}
    for group_key, n_units in group_counts.items():
        for sp_key, coeff in group_defs[group_key]["stoichiometric_ratio"].items():
            species_counts[sp_key] += int(n_units) * int(coeff)
    return species_counts

def get_atoms_per_molecule(species_cfg, species_key):
    atoms = species_cfg[species_key].get("atoms_per_molecule")
    if atoms is None:
        raise KeyError(
            f"Missing species.{species_key}.atoms_per_molecule required for system_sizing.mode='target_atoms'"
        )
    atoms = int(atoms)
    if atoms <= 0:
        raise ValueError(f"species.{species_key}.atoms_per_molecule must be > 0")
    return atoms

def estimate_box_size_nm(total_atoms, project_cfg, sizing_cfg):
    if sizing_cfg.get("estimate_box_from_atoms", False):
        atom_density = float(sizing_cfg.get("atom_number_density_atoms_per_nm3", 0.0))
        if atom_density <= 0:
            raise ValueError("system_sizing.atom_number_density_atoms_per_nm3 must be > 0 when estimate_box_from_atoms = true")
        volume_nm3 = float(total_atoms) / atom_density
        scaled_volume = project_cfg.get("box_scale_factor") * volume_nm3
        return (float(project_cfg.get("box_scale_factor",1.0)) * volume_nm3) ** (1.0 / 3.0), "estimated_and_scaled_from_atom_density"
    return float(project_cfg["box_size_nm"]), "project_settings.box_size_nm"

def resolve_system_size(cfg, species_cfg, order, group_defs, group_keys, component_ratio):
    project_cfg = cfg["project_settings"]
    sizing_cfg = cfg.get("system_sizing", {})
    mode = sizing_cfg.get("mode", "target_atoms")
    total_ratio_weight = sum(float(v) for v in component_ratio.values())
    if total_ratio_weight <= 0:
        raise ValueError("Component ratio must contain positive total weight")

    if mode == "reference_component":
        ref_component_key = sizing_cfg["reference_component_key"]
        ref_component_count = int(sizing_cfg["reference_component_count"])
        if ref_component_key not in group_keys:
            raise KeyError(
                f"system_sizing.reference_component_key '{ref_component_key}' is not a valid component. "
                f"Use one of: {', '.join(group_keys)}"
            )
        if ref_component_count < 0:
            raise ValueError("system_sizing.reference_component_count must be >= 0")
        ref_fraction = float(component_ratio[ref_component_key]) / total_ratio_weight
        if ref_fraction <= 0:
            raise ValueError(
                f"Component '{ref_component_key}' has zero weight in this screening.component_ratios entry"
            )
        n_total_components = int(round(ref_component_count / ref_fraction))
        component_counts = allocate_weighted_counts(n_total_components, component_ratio)
    elif mode == "target_atoms":
        target_atoms = int(sizing_cfg["target_atoms"])
        if target_atoms <= 0:
            raise ValueError("system_sizing.target_atoms must be > 0")
        avg_component_atoms = 0.0
        for group_key in group_keys:
            frac = float(component_ratio[group_key]) / total_ratio_weight
            group_atoms = sum(
                int(coeff) * get_atoms_per_molecule(species_cfg, sp_key)
                for sp_key, coeff in group_defs[group_key]["stoichiometric_ratio"].items()
            )
            avg_component_atoms += frac * group_atoms
        if avg_component_atoms <= 0:
            raise ValueError("Invalid sizing denominator. Check component ratios and species atom counts.")
        n_total_components = max(1, int(round(target_atoms / avg_component_atoms)))
        component_counts = allocate_weighted_counts(n_total_components, component_ratio)
    else:
        raise ValueError("system_sizing.mode must be one of: target_atoms, reference_component")
    group_counts = {k: int(component_counts[k]) for k in group_keys}
    counts = accumulate_species_counts_from_groups(group_counts, group_defs, order)

    total_atoms = None
    if mode == "target_atoms" or sizing_cfg.get("estimate_box_from_atoms", False):
        total_atoms = sum(counts[k] * get_atoms_per_molecule(species_cfg, k) for k in order)

    box_size_nm, box_source = (
        estimate_box_size_nm(total_atoms, project_cfg, sizing_cfg)
        if total_atoms is not None else
        (float(project_cfg["box_size_nm"]), "project_settings.box_size_nm")
    )

    return counts, {
        "mode": mode,
        "n_total_components": sum(group_counts.values()),
        "total_atoms": total_atoms,
        "box_size_nm": box_size_nm,
        "box_source": box_source,
        "component_counts": group_counts,
    }

def fmt_num(value):
    value = float(value)
    return str(int(value)) if value.is_integer() else f"{value:.4g}"

def build_case_label(temp, ratio_name, ff_name, component_ratio, active_itps, species_order, group_keys):
    comp_mix_label = "_".join(f"{k}{fmt_num(component_ratio[k])}" for k in group_keys)
    model_keys = species_order
    itp_label = "_".join(active_itps[k].split('.')[0] for k in model_keys)
    return f"T{fmt_num(temp)}_{ratio_name}_{ff_name}_{comp_mix_label}_{itp_label}"

def print_dry_run_case(label, temp, ratio_name, ff_name, counts, order, group_keys, component_ratio, active_itps, sizing_info):
    total_system = sum(counts[k] for k in order)
    component_weights = ", ".join(f"{k}={fmt_num(component_ratio[k])}" for k in group_keys)
    species_counts = ", ".join(f"{k}={counts[k]}" for k in order)
    print(f"{'='*100}\n[DRY-RUN] {label}\n{'='*100}")
    print(f"Temp={fmt_num(temp)} K\n Ratio={ratio_name}\n FF={ff_name}\n Component_ratio=({component_weights})")
    details = (
        f"Sizing_mode={sizing_info['mode']} -->  Total_groups={sizing_info['n_total_components']}, "
        f"box_nm={sizing_info['box_size_nm']:.3f} ({sizing_info['box_source']})"
    )
    print(details)
    print("  itp_map: " + ", ".join(f"{k}={active_itps[k]}" for k in order))
    component_counts = sizing_info.get("component_counts")
    if component_counts:
        print("  groups: " + ", ".join(f"{k}={component_counts[k]}" for k in group_keys))
    print(f"  counts: {species_counts}")
    totals_line = f"  totals: system={total_system}"
    if sizing_info["total_atoms"] is not None:
        totals_line += f", atoms={sizing_info['total_atoms']}"
    print(totals_line)

def append_file_to_aggregate_log(aggregate_path, source_path, section_title):
    if not os.path.exists(source_path):
        return
    try:
        with open(source_path, "r", encoding="utf-8", errors="ignore") as src:
            content = src.read().rstrip()
    except OSError:
        return
    with open(aggregate_path, "a", encoding="utf-8") as agg:
        agg.write(f"\n{'=' * 100}\n")
        agg.write(f"{section_title}\n")
        agg.write(f"Source: {source_path}\n")
        agg.write(f"{'=' * 100}\n")
        if content:
            agg.write(content)
            agg.write("\n")
        else:
            agg.write("[empty log]\n")

def rebuild_replica_grompp_log(rep_root):
    aggregate_path = os.path.join(rep_root, "grompp.log")
    os.makedirs(rep_root, exist_ok=True)
    with open(aggregate_path, "w", encoding="utf-8") as agg:
        agg.write("GROMACS aggregated replica setup/production log snapshot\n")
        agg.write("This file is regenerated on each launcher execution from existing stage logs,\n")
        agg.write("then appended by new setup/production grompp attempts for this replica.\n")

    stage_logs = [
        ("1_min", "insert-molecules.log"),
        ("1_min", "grompp_min.log"),
        ("2_press", "grompp_press.log"),
        ("3_anneal", "grompp_anneal.log"),
        ("4_prod", "grompp_prod.log"),
    ]
    replica_label = os.path.basename(rep_root)
    for stage, filename in stage_logs:
        source_path = os.path.join(rep_root, stage, filename)
        append_file_to_aggregate_log(
            aggregate_path,
            source_path,
            f"{replica_label} | {stage}/{filename}"
        )
    return aggregate_path

def build_case_contexts(master_combos, order, group_keys):
    contexts = []
    for ratio_entry, temp, ff_entry in master_combos:
        ratio_name = ratio_entry["name"]
        component_ratio = ratio_entry["weights"]
        ff_name = ff_entry["name"]
        active_itps = ff_entry["species_itps"]
        label = build_case_label(temp, ratio_name, ff_name, component_ratio, active_itps, order, group_keys)
        sys_path = os.path.abspath(f"data/{label}")
        contexts.append({
            "label": label,
            "sys_path": sys_path,
            "job_log_path": os.path.join(sys_path, "launch.log"),
            "ratio_entry": ratio_entry,
            "temp": temp,
            "ff_entry": ff_entry,
            "active_itps": active_itps,
        })
    return contexts

def prepare_replica_stage_inputs(rep_root, cfg, order, active_itps, sp, counts, temp, label):
    for stage in ["1_min", "2_press", "3_anneal", "4_prod"]:
        dest = os.path.join(rep_root, stage)
        os.makedirs(dest, exist_ok=True)
        for f in cfg["project_settings"]["global_ff_files"]:
            shutil.copy(f"inputs/{f}", f"{dest}/{f}")
        for k in order:
            shutil.copy(f"inputs/{active_itps[k]}", f"{dest}/{active_itps[k]}")
            shutil.copy(f"inputs/{sp[k]['gro']}", f"{dest}/{sp[k]['gro']}")

        with open(os.path.join(dest, "topol.top"), "w") as top:
            top.write('#include "forcefield.itp"\n')
            for k in order:
                top.write(f'#include "{active_itps[k]}"\n')
            top.write(f"\n[ system ]\n{label}\n\n[ molecules ]\n")
            for k in order:
                top.write(f"{sp[k]['resname']:<15} {counts[k]}\n")

        tdir, sim = cfg["project_settings"]["template_dir"], cfg["simulation_settings"]
        if stage == "1_min":
            shutil.copy(f"{tdir}/min.mdp", f"{dest}/min.mdp")
        elif stage == "2_press":
            modify_mdp(
                f"{tdir}/press.mdp",
                f"{dest}/press.mdp",
                {
                    "ref_t": sim["press"]["ref_t"],
                    "nsteps": sim["press"]["nsteps"],
                    "dt": sim["press"]["dt"],
                    "ref_p": sim["press"]["ref_p"],
                },
            )
        elif stage == "3_anneal":
            modify_mdp(
                f"{tdir}/anneal.mdp",
                f"{dest}/anneal.mdp",
                {
                    "ref_t": temp,
                    "ref_p": sim["anneal"]["ref_p"],
                    "nsteps": sim["anneal"]["nsteps"],
                    "dt": sim["anneal"]["dt"],
                    "annealing-time": sim["anneal"]["annealing_time"],
                    "annealing-temp": f"{sim['anneal']['temp_high']} {sim['anneal']['temp_high']} {temp} {temp}",
                },
            )
        elif stage == "4_prod":
            modify_mdp(
                f"{tdir}/nvt_run.mdp",
                f"{dest}/nvt_run.mdp",
                {
                    "ref_t": temp,
                    "nsteps": sim["prod"]["nsteps"],
                    "dt": sim["prod"]["dt"],
                },
            )

def build_initial_box_if_needed(rep_root, label, counts, order, sp, gmx_path, box, job_log, aggregate_log_path, replica_idx):
    min_dir = os.path.join(rep_root, "1_min")
    min_start_gro = os.path.join(min_dir, "start.gro")
    ins_log_abs = os.path.join(min_dir, "insert-molecules.log")
    curr = None
    if os.path.exists(min_start_gro):
        job_log(f"Replica {replica_idx}: found existing 1_min/start.gro at {min_start_gro}; skipping initial box rebuild.")
        return

    job_log(f"Replica {replica_idx}: building initial box in 1_min (target box={box:.3f} nm).")
    with open(ins_log_abs, "w") as l:
        l.write("=== BOX LOG ===\n")

    for k in order:
        if counts[k] <= 0:
            with open(ins_log_abs, "a") as l:
                l.write(f"\n--- Skipping {k} (nmol=0) ---\n")
            job_log(f"Replica {replica_idx}: skipped species '{k}' during insert-molecules (nmol=0).")
            continue
        out = f"tmp_{k}_{replica_idx}.gro"
        cmd = (
            f"mpirun -np 1 {gmx_path} insert-molecules "
            f"{'-box '+str(box)+' '+str(box)+' '+str(box) if curr is None else '-f '+curr} "
            f"-ci inputs/{sp[k]['gro']} -nmol {counts[k]} -o {out} -seed {replica_idx*123}"
        )
        job_log(f"Replica {replica_idx}: insert-molecules for '{k}' (nmol={counts[k]}).")
        res = subprocess.run(cmd.split(), capture_output=True, text=True)
        with open(ins_log_abs, "a") as l:
            l.write(f"\n--- Adding {k} ---\n" + res.stdout + res.stderr)
        if not os.path.exists(out):
            print(f"FAILED build for {label}. Check {ins_log_abs}")
            sys.exit(1)
        curr = out

    if curr is None:
        print(f"FAILED build for {label}: no molecules were inserted. Check {ins_log_abs}")
        job_log(f"Replica {replica_idx}: build failed, no molecules inserted. Check {ins_log_abs}")
        sys.exit(1)

    shutil.move(curr, min_start_gro)
    job_log(f"Replica {replica_idx}: built start.gro at {min_start_gro}")
    append_file_to_aggregate_log(
        aggregate_log_path,
        ins_log_abs,
        f"Replica {replica_idx} | 1_min/insert-molecules.log"
    )
    for f in os.listdir("."):
        if f.startswith("tmp_") and f.endswith(".gro"):
            os.remove(f)

def replica_has_existing_state(rep_root):
    if not os.path.exists(rep_root):
        return False
    try:
        return any(True for _ in os.scandir(rep_root))
    except OSError:
        return True

def write_setup_sh(path, cfg, abs_rep_path, aggregate_log_path):
    p, s = cfg["project_settings"], cfg["slurm_settings"]
    module_block = "\n".join(p.get("gmx_modules", []))
    analysis_script = os.path.join(os.getcwd(), p['scripts_dir'], p['density_analysis_script'])
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=S_{p['system_name']}
#SBATCH --account={s['account']}
#SBATCH --partition={s['partition_setup']}
#SBATCH --qos={s['qos_setup']}
#SBATCH --time={s['time_setup']}
#SBATCH --nodes={s.get('nodes_setup', 1)}
#SBATCH --ntasks-per-node={s.get('ntasks_per_node_setup', 56)}
#SBATCH --cpus-per-task={s.get('cpus_per_task_setup', 2)}
#SBATCH --output={abs_rep_path}/setup.out
#SBATCH --error={abs_rep_path}/setup.err

set -euo pipefail

{module_block}
source ~/miniconda3/bin/activate {p['conda_env']}

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

BASE="{abs_rep_path}"
GMX="{p.get('gmx_executable_path', 'gmx_mpi')}"
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

# 1. MIN
cd $BASE/1_min
if [ ! -f "min_out.gro" ]; then
    if ! srun -n 1 $GMX grompp -f min.mdp -c start.gro -p topol.top -o min.tpr > $BASE/1_min/grompp_min.log 2>&1; then
        append_stage_log "Replica $(basename "$BASE") | 1_min/grompp_min.log" "$BASE/1_min/grompp_min.log"
        echo "ERROR: grompp failed for minimization. See $BASE/1_min/grompp_min.log"
        exit 1
    fi
    append_stage_log "Replica $(basename "$BASE") | 1_min/grompp_min.log" "$BASE/1_min/grompp_min.log"
    srun $GMX mdrun -deffnm min -c min_out.gro
else
    echo "Minimization already completed. Skipping."
fi

# 2. PRESS
cd $BASE/2_press
if [ ! -f "press_out.gro" ]; then
    if ! srun -n 1 $GMX grompp -f press.mdp -c ../1_min/min_out.gro -p topol.top -o press.tpr > $BASE/2_press/grompp_press.log 2>&1; then
        append_stage_log "Replica $(basename "$BASE") | 2_press/grompp_press.log" "$BASE/2_press/grompp_press.log"
        echo "ERROR: grompp failed for pressure equilibration. See $BASE/2_press/grompp_press.log"
        exit 1
    fi
    append_stage_log "Replica $(basename "$BASE") | 2_press/grompp_press.log" "$BASE/2_press/grompp_press.log"
    srun $GMX mdrun -deffnm press -c press_out.gro -rdd 1.2 -pin on
else
    echo "Pressure equilibration already completed. Skipping."
fi

# 3. ANNEAL
cd $BASE/3_anneal
if [ ! -f "anneal_out.gro" ]; then
    if ! srun -n 1 $GMX grompp -f anneal.mdp -c ../2_press/press_out.gro -p topol.top -o anneal.tpr > $BASE/3_anneal/grompp_anneal.log 2>&1; then
        append_stage_log "Replica $(basename "$BASE") | 3_anneal/grompp_anneal.log" "$BASE/3_anneal/grompp_anneal.log"
        echo "ERROR: grompp failed for annealing. See $BASE/3_anneal/grompp_anneal.log"
        exit 1
    fi
    append_stage_log "Replica $(basename "$BASE") | 3_anneal/grompp_anneal.log" "$BASE/3_anneal/grompp_anneal.log"
    srun $GMX mdrun -deffnm anneal -c anneal_out.gro -pin on
else
    echo "Annealing already completed. Skipping."
fi

# 4. ANALYSIS
cd $BASE/3_anneal
if [ ! -f "$BASE/4_prod/start.gro" ]; then
    python3 {analysis_script} anneal.edr --num_replicas 1 --time_crop {p['density_analysis_time_crop']}
    if [ -f "start_replica_1.gro" ]; then
        mv start_replica_1.gro $BASE/4_prod/start.gro
    else
        echo "ERROR: density analysis did not produce start_replica_1.gro"
        exit 1
    fi
else
    echo "Analysis and start.gro for production already exist. Skipping."
fi
"""
    with open(os.path.join(path, "run_setup.sh"), "w") as f: f.write(script_content)

def write_prod_sh(path, cfg, abs_rep_path, chunk_idx, aggregate_log_path):
    p, s = cfg["project_settings"], cfg["slurm_settings"]
    module_block = "\n".join(p.get("gmx_modules", []))
    prod_path = os.path.join(abs_rep_path, "4_prod")
    gmx_exe = p.get('gmx_executable_path', 'gmx_mpi')
    script_content = f"""#!/bin/bash
#SBATCH --job-name=P{chunk_idx}_{p['system_name']}
#SBATCH --account={s['account']}
#SBATCH --partition={s['partition_prod']}
#SBATCH --qos={s['qos_prod']}
#SBATCH --time={s['time_prod']}
#SBATCH --nodes={s.get('nodes_prod', 1)}
#SBATCH --ntasks-per-node={s.get('ntasks_per_node_prod', 56)}
#SBATCH --cpus-per-task={s.get('cpus_per_task_prod', 2)}
#SBATCH --output={prod_path}/chunk_{chunk_idx}.out
#SBATCH --error={prod_path}/chunk_{chunk_idx}.err

set -euo pipefail

{module_block}

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

GMX="{p.get('gmx_executable_path', 'gmx_mpi')}"
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

cd {prod_path}
# Select a TPR file dynamically (supports names like prun.tpr, production.tpr, etc.).
TPR_FILE=""
if [ -f "nvt_run.tpr" ]; then
    TPR_FILE="nvt_run.tpr"
else
    shopt -s nullglob
    TPR_CANDIDATES=( *.tpr )
    shopt -u nullglob
    if [ "${{#TPR_CANDIDATES[@]}}" -gt 0 ]; then
        TPR_FILE="${{TPR_CANDIDATES[0]}}"
    fi
fi

if [ -n "$TPR_FILE" ]; then
    DEFFNM="${{TPR_FILE%.tpr}}"
    CPT_FILE="${{DEFFNM}}.cpt"
else
    DEFFNM="nvt_run"
    CPT_FILE="nvt_run.cpt"
fi

if [ -f "$CPT_FILE" ]; then
    if [ ! -f "$TPR_FILE" ]; then
        echo "ERROR: $CPT_FILE exists but no matching .tpr file was found. Cannot restart production."
        exit 1
    fi
    srun $GMX mdrun -v -deffnm "$DEFFNM" -s "$TPR_FILE" -cpi "$CPT_FILE" -append -maxh 71.7 -pin on
else
    if [ -n "$TPR_FILE" ]; then
        srun $GMX mdrun -v -deffnm "$DEFFNM" -s "$TPR_FILE" -maxh 71.7 -pin on
    else
        if ! srun -n 1 $GMX grompp -f nvt_run.mdp -c start.gro -p topol.top -o nvt_run.tpr > {prod_path}/grompp_prod.log 2>&1; then
            append_stage_log "Replica $(basename "$(dirname "$PWD")") | 4_prod/grompp_prod.log" "{prod_path}/grompp_prod.log"
            echo "ERROR: grompp failed for production. See {prod_path}/grompp_prod.log"
            exit 1
        fi
        append_stage_log "Replica $(basename "$(dirname "$PWD")") | 4_prod/grompp_prod.log" "{prod_path}/grompp_prod.log"
        srun $GMX mdrun -v -deffnm nvt_run -s nvt_run.tpr -maxh 71.7 -pin on
    fi
fi
"""
    with open(os.path.join(path, f"run_prod_{chunk_idx}.sh"), "w") as f: f.write(script_content)

def make_job_logger(sys_path, label):
    os.makedirs(sys_path, exist_ok=True)
    log_path = os.path.join(sys_path, "launch.log")

    def _log(msg):
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(log_path, "a", encoding="utf-8") as f:
            f.write("\n")
            for line in str(msg).splitlines():
                f.write(f"[{ts}] {line}\n")

    _log("=" * 100)
    _log(f"LAUNCH WORKFLOW START: {label}")
    _log("=" * 100)
    return _log, log_path

def format_final_molar_ratio(component_counts, group_keys):
    values = [int(component_counts.get(k, 0)) for k in group_keys]
    positive = [v for v in values if v > 0]
    if not positive:
        return "all components are zero"
    gcd_val = positive[0]
    for v in positive[1:]:
        gcd_val = math.gcd(gcd_val, v)
    reduced = {k: (int(component_counts.get(k, 0)) // gcd_val) for k in group_keys}
    total = sum(values)
    fractions = {k: (int(component_counts.get(k, 0)) / total if total > 0 else 0.0) for k in group_keys}
    reduced_str = ":".join(str(reduced[k]) for k in group_keys)
    frac_str = ", ".join(f"{k}={fractions[k]:.6f}" for k in group_keys)
    return f"reduced={reduced_str} ({', '.join(group_keys)}), mole_fractions=({frac_str})"

def is_setup_complete(rep_root):
    required = [
        os.path.join(rep_root, "1_min", "min_out.gro"),
        os.path.join(rep_root, "2_press", "press_out.gro"),
        os.path.join(rep_root, "3_anneal", "anneal_out.gro"),
        os.path.join(rep_root, "4_prod", "start.gro"),
    ]
    return all(os.path.exists(p) for p in required)

def is_production_complete(prod_dir):
    log_path = os.path.join(prod_dir, "nvt_run.log")
    if not os.path.exists(log_path):
        return False
    try:
        with open(log_path, "r", errors="ignore") as f:
            content = f.read()
    except OSError:
        return False
    return "Finished mdrun" in content

def get_target_prod_steps(prod_dir, cfg):
    mdp_path = os.path.join(prod_dir, "nvt_run.mdp")
    if os.path.exists(mdp_path):
        try:
            with open(mdp_path, "r", errors="ignore") as f:
                for line in f:
                    cleaned = line.split(";")[0].strip()
                    if "=" not in cleaned:
                        continue
                    key, val = [x.strip() for x in cleaned.split("=", 1)]
                    if key == "nsteps":
                        return int(float(val))
        except (OSError, ValueError):
            pass
    return int(cfg["simulation_settings"]["prod"]["nsteps"])

def parse_chunk_err_status(prod_dir, target_steps):
    """
    Inspect chunk_*.err files in sequence (chunk_1, chunk_2, ...)
    and estimate whether production reached the requested total simulation steps.
    """
    completed_chunks = 0
    max_step = 0
    reached_target = False
    step_pattern = re.compile(r"\bstep\s+([0-9]+)\b")

    idx = 1
    last_seen_chunk_idx = 0
    while True:
        err_path = os.path.join(prod_dir, f"chunk_{idx}.err")
        if not os.path.exists(err_path):
            break
        last_seen_chunk_idx = idx
        try:
            with open(err_path, "r", errors="ignore") as f:
                content = f.read()
        except OSError:
            break

        steps = [int(m.group(1)) for m in step_pattern.finditer(content)]
        if steps:
            max_step = max(max_step, max(steps))
        if max_step >= target_steps:
            reached_target = True
            completed_chunks += 1
            break

        # A finished chunk that did not reach target means continuation was needed.
        if ("Writing final coordinates." in content) or ("Finished mdrun" in content):
            completed_chunks += 1
            idx += 1
            continue

        # Chunk exists but did not finish: stop scanning.
        break

    return {
        "completed_chunks": completed_chunks,
        "max_step": max_step,
        "reached_target": reached_target,
        "next_chunk_idx": max(1, completed_chunks + 1),
        "last_seen_chunk_idx": last_seen_chunk_idx,
    }

def get_next_unused_chunk_idx(prod_dir):
    pattern = re.compile(r"^chunk_(\d+)\.(?:out|err)$")
    max_idx = 0
    try:
        for name in os.listdir(prod_dir):
            match = pattern.match(name)
            if match:
                max_idx = max(max_idx, int(match.group(1)))
    except OSError:
        return 1
    return max_idx + 1 if max_idx > 0 else 1

def normalize_slurm_state(state):
    if not state:
        return "unknown"
    upper = state.strip().upper()
    if upper in {"PENDING", "CONFIGURING", "RESV_DEL_HOLD", "REQUEUE_HOLD", "REQUEUED", "SPECIAL_EXIT"}:
        return "queued"
    if upper in {"RUNNING", "COMPLETING", "STAGE_OUT", "SUSPENDED"}:
        return "running"
    if upper in {"COMPLETED"}:
        return "completed"
    if upper in {"FAILED", "NODE_FAIL", "BOOT_FAIL", "OUT_OF_MEMORY"}:
        return "failed"
    if upper in {"CANCELLED", "DEADLINE", "PREEMPTED", "REVOKED"} or upper.startswith("CANCELLED"):
        return "cancelled"
    if upper in {"TIMEOUT"}:
        return "timeout"
    return upper.lower()

def parse_replica_job_ids(log_path, replica_idx):
    result = {"setup_job_id": None, "prod_job_id": None}
    if not os.path.exists(log_path):
        return result
    setup_pattern = re.compile(
        rf"Replica {replica_idx}: submitted setup job with id ([^\s]+)"
    )
    prod_pattern = re.compile(
        rf"Replica {replica_idx}: submitted production chunk \d+ with job id ([^\s]+)"
    )
    try:
        with open(log_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                setup_match = setup_pattern.search(line)
                if setup_match:
                    result["setup_job_id"] = setup_match.group(1)
                prod_match = prod_pattern.search(line)
                if prod_match:
                    result["prod_job_id"] = prod_match.group(1)
    except OSError:
        return result
    return result

def query_squeue_status(job_ids):
    if not job_ids or not shutil.which("squeue"):
        return {}
    try:
        res = subprocess.run(
            [
                "squeue",
                "--noheader",
                "--format=%i|%T",
                "--jobs",
                ",".join(job_ids),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return {}
    if res.returncode != 0:
        return {}
    statuses = {}
    for line in res.stdout.splitlines():
        parts = line.strip().split("|", 1)
        if len(parts) != 2:
            continue
        job_id, state = parts
        statuses[job_id.strip()] = normalize_slurm_state(state)
    return statuses

def query_sacct_status(job_ids):
    if not job_ids or not shutil.which("sacct"):
        return {}
    try:
        res = subprocess.run(
            [
                "sacct",
                "-n",
                "-P",
                "-X",
                "-j",
                ",".join(job_ids),
                "--format=JobIDRaw,State",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return {}
    if res.returncode != 0:
        return {}
    statuses = {}
    for line in res.stdout.splitlines():
        parts = line.strip().split("|")
        if len(parts) < 2:
            continue
        job_id, state = parts[0].strip(), parts[1].strip()
        if not job_id or "." in job_id:
            continue
        statuses[job_id] = normalize_slurm_state(state)
    return statuses

def collect_slurm_status_map(job_ids):
    unique_ids = [job_id for job_id in dict.fromkeys(job_ids) if job_id]
    squeue_status = query_squeue_status(unique_ids)
    sacct_needed = [job_id for job_id in unique_ids if job_id not in squeue_status]
    sacct_status = query_sacct_status(sacct_needed)
    status_map = {}
    for job_id in unique_ids:
        if job_id in squeue_status:
            status_map[job_id] = squeue_status[job_id]
        elif job_id in sacct_status:
            status_map[job_id] = sacct_status[job_id]
    return status_map

def status_flag(condition):
    return "yes" if condition else "no"

def build_progress_row(case_label, rep_idx, rep_root, cfg, launch_log_path, slurm_status_map):
    min_dir = os.path.join(rep_root, "1_min")
    press_dir = os.path.join(rep_root, "2_press")
    anneal_dir = os.path.join(rep_root, "3_anneal")
    prod_dir = os.path.join(rep_root, "4_prod")

    target_steps = get_target_prod_steps(prod_dir, cfg)
    chunk_status = parse_chunk_err_status(prod_dir, target_steps)
    job_ids = parse_replica_job_ids(launch_log_path, rep_idx)
    tracked_job_id = job_ids["prod_job_id"] or job_ids["setup_job_id"]
    slurm_status = slurm_status_map.get(tracked_job_id)
    if not slurm_status:
        if chunk_status["reached_target"]:
            slurm_status = "completed"
        elif tracked_job_id:
            slurm_status = "unknown"
        else:
            slurm_status = "not_submitted"

    row = {
        "simulation": case_label,
        "replica": rep_idx,
        "slurm_job_id": tracked_job_id or "",
        "slurm_status": slurm_status,
        "insert_molecules": status_flag(os.path.exists(os.path.join(min_dir, "start.gro"))),
        "grompp_min": status_flag(os.path.exists(os.path.join(min_dir, "min.tpr"))),
        "min_run": status_flag(os.path.exists(os.path.join(min_dir, "min_out.gro"))),
        "grompp_press": status_flag(os.path.exists(os.path.join(press_dir, "press.tpr"))),
        "press_run": status_flag(os.path.exists(os.path.join(press_dir, "press_out.gro"))),
        "grompp_anneal": status_flag(os.path.exists(os.path.join(anneal_dir, "anneal.tpr"))),
        "anneal_run": status_flag(os.path.exists(os.path.join(anneal_dir, "anneal_out.gro"))),
        "analysis_start_gro": status_flag(os.path.exists(os.path.join(prod_dir, "start.gro"))),
        "grompp_prod": status_flag(os.path.exists(os.path.join(prod_dir, "nvt_run.tpr"))),
        "prod_checkpoint": status_flag(os.path.exists(os.path.join(prod_dir, "nvt_run.cpt"))),
        "prod_target_reached": status_flag(chunk_status["reached_target"]),
        "max_step_seen": chunk_status["max_step"],
        "target_steps": target_steps,
        "next_chunk_idx": max(
            chunk_status["next_chunk_idx"],
            get_next_unused_chunk_idx(prod_dir)
        ),
    }
    return row

def write_progress_table(rows, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, "simulation_progress.csv")
    md_path = os.path.join(output_dir, "simulation_progress.md")
    fieldnames = [
        "simulation",
        "replica",
        "slurm_job_id",
        "slurm_status",
        "insert_molecules",
        "grompp_min",
        "min_run",
        "grompp_press",
        "press_run",
        "grompp_anneal",
        "anneal_run",
        "analysis_start_gro",
        "grompp_prod",
        "prod_checkpoint",
        "prod_target_reached",
        "max_step_seen",
        "target_steps",
        "next_chunk_idx",
    ]

    with open(csv_path, "w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    headers = fieldnames
    header_labels = {
        "simulation": "simulation",
        "replica": "replica",
        "slurm_job_id": "slurm job id",
        "slurm_status": "slurm status",
        "insert_molecules": "insert-molecules",
        "grompp_min": "grompp min",
        "min_run": "min run",
        "grompp_press": "grompp press",
        "press_run": "press run",
        "grompp_anneal": "grompp anneal",
        "anneal_run": "anneal run",
        "analysis_start_gro": "analysis/start.gro",
        "grompp_prod": "grompp prod",
        "prod_checkpoint": "prod cpt",
        "prod_target_reached": "target steps reached",
        "max_step_seen": "max step seen",
        "target_steps": "target steps",
        "next_chunk_idx": "next chunk",
    }
    with open(md_path, "w", encoding="utf-8") as md_file:
        md_file.write("# Simulation Progress\n\n")
        md_file.write("Automatically generated by `launch_gromacs.py`.\n\n")
        md_file.write("| " + " | ".join(header_labels[h] for h in headers) + " |\n")
        md_file.write("| " + " | ".join("---" for _ in headers) + " |\n")
        for row in rows:
            md_file.write("| " + " | ".join(str(row[h]) for h in headers) + " |\n")

    return csv_path, md_path

def collect_progress_snapshot(case_contexts, cfg):
    slurm_job_ids = []
    progress_row_refs = []
    for ctx in case_contexts:
        for replica_idx in range(1, cfg['project_settings']['num_replicas'] + 1):
            rep_root = os.path.join(ctx["sys_path"], f"rep_{replica_idx}")
            progress_row_refs.append((ctx["label"], replica_idx, rep_root))
            job_ids = parse_replica_job_ids(ctx["job_log_path"], replica_idx)
            if job_ids["setup_job_id"]:
                slurm_job_ids.append(job_ids["setup_job_id"])
            if job_ids["prod_job_id"]:
                slurm_job_ids.append(job_ids["prod_job_id"])

    slurm_status_map = collect_slurm_status_map(slurm_job_ids)
    log_path_by_label = {ctx["label"]: ctx["job_log_path"] for ctx in case_contexts}
    return [
        build_progress_row(
            label,
            replica_idx,
            rep_root,
            cfg,
            log_path_by_label[label],
            slurm_status_map,
        )
        for label, replica_idx, rep_root in progress_row_refs
    ]

def submit_sbatch(script_path, dependency_job_id=None):
    cmd = ["sbatch", "--parsable"]
    if dependency_job_id:
        cmd.append(f"--dependency=afterany:{dependency_job_id}")
    cmd.append(script_path)
    try:
        res = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise RuntimeError(f"Failed to execute sbatch for '{script_path}': {exc}") from exc

    stdout = (res.stdout or "").strip()
    stderr = (res.stderr or "").strip()
    if res.returncode != 0:
        raise RuntimeError(
            f"sbatch failed for '{script_path}' with return code {res.returncode}. "
            f"stderr: {stderr or '[empty]'}"
        )
    if not stdout:
        raise RuntimeError(
            f"sbatch returned no job id for '{script_path}'. "
            f"stderr: {stderr or '[empty]'}"
        )
    return stdout

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", nargs='?', default="config.toml")
    parser.add_argument("--dry-run", action="store_true", help="Print computed species counts and exit without running GROMACS or submitting jobs")
    parser.add_argument(
        "--track-progress",
        action="store_true",
        help="Do not submit or modify simulations; only refresh the progress tables from existing files/logs and Slurm state",
    )
    args = parser.parse_args()

    if args.dry_run and args.track_progress:
        print("ERROR: --dry-run and --track-progress cannot be used together.")
        sys.exit(1)

    if not os.path.exists(args.config): create_default_config(args.config)
    with open(args.config, "rb") as f: cfg = tomllib.load(f)

    scr = cfg["screening"]
    sp = normalize_species_config(cfg["species"])
    def to_list(val): return val if isinstance(val, list) else [val]
    group_defs, group_keys, active_species_keys = discover_groups(cfg, sp)
    order = active_species_keys

    # --- 1. DEFINE SCREENING AXES ---
    component_ratios = normalize_component_ratios(scr, group_keys)
    force_field_sets = normalize_force_field_sets(scr, order)

    # --- 2. CROSS COMBINATIONS ---
    # Cross (Component Ratio) x (Temperature) x (Force-Field Set)
    master_combos = list(itertools.product(component_ratios, to_list(scr["target_temps"]), force_field_sets))

    print(f"Setting {len(master_combos)} total variations...")

    case_contexts = build_case_contexts(master_combos, order, group_keys)

    if args.track_progress:
        progress_rows = collect_progress_snapshot(case_contexts, cfg)
        progress_csv_path, progress_md_path = write_progress_table(
            progress_rows,
            os.path.abspath("data")
        )
        print(f"Track-progress mode: refreshed {progress_md_path} and {progress_csv_path}")
        return

    if not args.dry_run:
        # --- GMX DETECTION ---
        gmx_path = cfg["project_settings"].get("gmx_executable_path")
        if not gmx_path or not os.path.exists(gmx_path): gmx_path = shutil.which("gmx_mpi")
        if not gmx_path: print("ERROR: GMX not found"); sys.exit(1)
        

    progress_row_refs = []
    slurm_job_ids = []

    for ctx in case_contexts:
        ratio_entry = ctx["ratio_entry"]
        temp = ctx["temp"]
        ff_entry = ctx["ff_entry"]
        active_itps = ctx["active_itps"]
        ratio_name = ratio_entry["name"]
        ff_name = ff_entry["name"]
        label = ctx["label"]
        sys_path = ctx["sys_path"]
        job_log, job_log_path = make_job_logger(sys_path, label)
        job_log(f"Case details: ratio={ratio_name}, ff={ff_name}, temp={fmt_num(temp)} K")
        job_log(f"Log file path: {job_log_path}")
        
        counts, sizing_info = resolve_system_size(
            cfg, sp, order, group_defs, group_keys,
            ratio_entry["weights"]
        )
        job_log(
            "Resolved system size: "
            f"mode={sizing_info['mode']}, total_groups={sizing_info['n_total_components']}, "
            f"box_nm={sizing_info['box_size_nm']:.3f}, "
            f"species_counts={', '.join(f'{k}={counts[k]}' for k in order)}"
        )
        component_counts = sizing_info.get("component_counts", {})
        job_log(
            "Final achieved molar ratio from rounded group counts: "
            + format_final_molar_ratio(component_counts, group_keys)
        )
        for r in range(1, cfg['project_settings']['num_replicas'] + 1):
            job_ids = parse_replica_job_ids(job_log_path, r)
            if job_ids["setup_job_id"]:
                slurm_job_ids.append(job_ids["setup_job_id"])
            if job_ids["prod_job_id"]:
                slurm_job_ids.append(job_ids["prod_job_id"])
        ctx["job_log_path"] = job_log_path

        if args.dry_run:
            job_log("Dry-run enabled. No file generation, no GROMACS, no Slurm submission.")
            print_dry_run_case(
                label,
                temp,
                ratio_name,
                ff_name,
                counts,
                order,
                group_keys,
                ratio_entry["weights"],
                active_itps,
                sizing_info,
            )
            continue

        for r in range(1, cfg['project_settings']['num_replicas'] + 1):
            rep_root = os.path.join(sys_path, f"rep_{r}")
            if replica_has_existing_state(rep_root):
                msg = (
                    f"Replica {r}: existing replica folder detected at {rep_root}. "
                    "Initial-launch mode will not modify or resubmit existing replicas."
                )
                print(f"{label} R{r}: existing replica detected. Skipping initial submission.")
                job_log(msg)
                progress_row_refs.append((label, r, rep_root))
                continue
            grompp_log_path = rebuild_replica_grompp_log(rep_root)
            job_log("-" * 70)
            job_log(f"Replica {r}")
            job_log("-" * 70)
            job_log(f"Replica {r}: preparing stage folders and topology/MDP inputs in {rep_root}")
            job_log(f"Replica {r}: aggregated grompp log path: {grompp_log_path}")
            prepare_replica_stage_inputs(rep_root, cfg, order, active_itps, sp, counts, temp, label)
            build_initial_box_if_needed(
                rep_root,
                label,
                counts,
                order,
                sp,
                gmx_path,
                sizing_info["box_size_nm"],
                job_log,
                grompp_log_path,
                r,
            )

            # --- PRE-SUBMISSION CHECKS ---
            prod_dir = os.path.join(rep_root, "4_prod")
            setup_done = is_setup_complete(rep_root)
            has_start = os.path.exists(os.path.join(prod_dir, "start.gro"))
            has_checkpoint = os.path.exists(os.path.join(prod_dir, "nvt_run.cpt"))
            target_steps = get_target_prod_steps(prod_dir, cfg)
            chunk_status = parse_chunk_err_status(prod_dir, target_steps)
            job_log("Pre-submission checks")
            job_log("-" * 70)
            job_log(
                f"Replica {r}: pre-submit checks -> setup_done={setup_done}, has_start={has_start}, "
                f"has_checkpoint={has_checkpoint}, target_steps={target_steps}, "
                f"max_step_seen={chunk_status['max_step']}, completed_chunks={chunk_status['completed_chunks']}, "
                f"next_chunk_idx={chunk_status['next_chunk_idx']}, "
                f"last_seen_chunk_idx={chunk_status['last_seen_chunk_idx']}"
            )

            if chunk_status["reached_target"]:
                msg = (
                    f"Skipping {label} R{r}: production reached target steps "
                    f"({chunk_status['max_step']} / {target_steps})."
                )
                print(msg)
                job_log(f"Replica {r}: {msg}")
                progress_row_refs.append((label, r, rep_root))
                continue

            prev = None
            if not setup_done:
                write_setup_sh(rep_root, cfg, rep_root, grompp_log_path)
                try:
                    sid = submit_sbatch(os.path.join(rep_root, "run_setup.sh"))
                except RuntimeError as exc:
                    job_log(f"Replica {r}: ERROR submitting setup job: {exc}")
                    print(f"ERROR submitting setup job for {label} R{r}: {exc}")
                    sys.exit(1)
                if sid:
                    slurm_job_ids.append(sid)
                print(f"Launched {label} R{r}: Setup {sid}")
                job_log(f"Replica {r}: submitted setup job with id {sid}")
                prev = sid
            else:
                print(f"{label} R{r}: setup already complete. Skipping setup submission.")
                job_log(f"Replica {r}: setup already complete, setup submission skipped.")

            # Production launch policy:
            # Keep support for multi-chunk submission per launcher execution.
            # The run script reuses existing nvt_run.tpr and never relies on
            # backup-renamed files.
            if has_start and not has_checkpoint:
                start_chunk_idx = get_next_unused_chunk_idx(prod_dir)
                prod_jobs = cfg["project_settings"]["num_prod_chunks"]
            elif has_start and has_checkpoint:
                start_chunk_idx = max(
                    chunk_status["next_chunk_idx"],
                    get_next_unused_chunk_idx(prod_dir)
                )
                prod_jobs = 1
            else:
                start_chunk_idx = get_next_unused_chunk_idx(prod_dir)
                prod_jobs = cfg["project_settings"]["num_prod_chunks"]
            job_log(
                f"Replica {r}: production policy -> start_chunk_idx={start_chunk_idx}, "
                f"prod_jobs_to_submit={prod_jobs}"
            )

            for c in range(start_chunk_idx, start_chunk_idx + prod_jobs):
                write_prod_sh(rep_root, cfg, rep_root, c, grompp_log_path)
                try:
                    prev = submit_sbatch(
                        os.path.join(rep_root, f"run_prod_{c}.sh"),
                        dependency_job_id=prev,
                    )
                except RuntimeError as exc:
                    job_log(f"Replica {r}: ERROR submitting production chunk {c}: {exc}")
                    print(f"ERROR submitting production chunk {c} for {label} R{r}: {exc}")
                    sys.exit(1)
                if prev:
                    slurm_job_ids.append(prev)
                job_log(f"Replica {r}: submitted production chunk {c} with job id {prev}")
            print(
                f"Launched {label} R{r}: {prod_jobs} production job(s) submitted "
                f"(starting at chunk {start_chunk_idx})."
            )
            job_log(
                f"Replica {r}: production submissions complete "
                f"(count={prod_jobs}, starting_chunk={start_chunk_idx})"
            )

            progress_row_refs.append((label, r, rep_root))

    if not progress_row_refs:
        for ctx in case_contexts:
            for r in range(1, cfg['project_settings']['num_replicas'] + 1):
                rep_root = os.path.join(ctx["sys_path"], f"rep_{r}")
                progress_row_refs.append((ctx["label"], r, rep_root))

    slurm_status_map = collect_slurm_status_map(slurm_job_ids)
    log_path_by_label = {ctx["label"]: ctx["job_log_path"] for ctx in case_contexts}
    progress_rows = [
        build_progress_row(
            label,
            replica_idx,
            rep_root,
            cfg,
            log_path_by_label[label],
            slurm_status_map,
        )
        for label, replica_idx, rep_root in progress_row_refs
    ]

    progress_csv_path, progress_md_path = write_progress_table(
        progress_rows,
        os.path.abspath("data")
    )
    print(f"Simulation progress table written to {progress_md_path} and {progress_csv_path}")

if __name__ == "__main__": main()

