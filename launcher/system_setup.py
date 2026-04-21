import os
import shutil
import subprocess
import sys

from .mdp import modify_mdp


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


def prepare_replica_stage_inputs(rep_root, cfg, order, active_itps, species_cfg, counts, temp, label):
    topology_include = cfg["project_settings"].get("topology_forcefield_include", "forcefield.itp")
    for stage in ["1_min", "2_press", "3_anneal", "4_prod"]:
        dest = os.path.join(rep_root, stage)
        os.makedirs(dest, exist_ok=True)
        for filename in cfg["project_settings"]["global_ff_files"]:
            shutil.copy(f"inputs/{filename}", f"{dest}/{filename}")
        for key in order:
            shutil.copy(f"inputs/{active_itps[key]}", f"{dest}/{active_itps[key]}")
            shutil.copy(f"inputs/{species_cfg[key]['gro']}", f"{dest}/{species_cfg[key]['gro']}")

        with open(os.path.join(dest, "topol.top"), "w", encoding="utf-8") as top:
            top.write(f'#include "{topology_include}"\n')
            for key in order:
                top.write(f'#include "{active_itps[key]}"\n')
            top.write(f"\n[ system ]\n{label}\n\n[ molecules ]\n")
            for key in order:
                top.write(f"{species_cfg[key]['resname']:<15} {counts[key]}\n")

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


def build_initial_box_if_needed(rep_root, label, counts, order, species_cfg, gmx_path, box, job_log, aggregate_log_path, replica_idx):
    min_dir = os.path.join(rep_root, "1_min")
    min_start_gro = os.path.join(min_dir, "start.gro")
    ins_log_abs = os.path.join(min_dir, "insert-molecules.log")
    curr = None
    if os.path.exists(min_start_gro):
        job_log(
            f"Replica {replica_idx}: found existing 1_min/start.gro at {min_start_gro}; skipping initial box rebuild."
        )
        return

    job_log(f"Replica {replica_idx}: building initial box in 1_min (target box={box:.3f} nm).")
    with open(ins_log_abs, "w", encoding="utf-8") as log_file:
        log_file.write("=== BOX LOG ===\n")

    for key in order:
        if counts[key] <= 0:
            with open(ins_log_abs, "a", encoding="utf-8") as log_file:
                log_file.write(f"\n--- Skipping {key} (nmol=0) ---\n")
            job_log(f"Replica {replica_idx}: skipped species '{key}' during insert-molecules (nmol=0).")
            continue
        out = f"tmp_{key}_{replica_idx}.gro"
        cmd = (
            f"mpirun -np 1 {gmx_path} insert-molecules "
            f"{'-box ' + str(box) + ' ' + str(box) + ' ' + str(box) if curr is None else '-f ' + curr} "
            f"-ci inputs/{species_cfg[key]['gro']} -nmol {counts[key]} -o {out} -seed {replica_idx * 123}"
        )
        job_log(f"Replica {replica_idx}: insert-molecules for '{key}' (nmol={counts[key]}).")
        res = subprocess.run(cmd.split(), capture_output=True, text=True)
        with open(ins_log_abs, "a", encoding="utf-8") as log_file:
            log_file.write(f"\n--- Adding {key} ---\n{res.stdout}{res.stderr}")
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
    for filename in os.listdir("."):
        if filename.startswith("tmp_") and filename.endswith(".gro"):
            os.remove(filename)


def replica_has_existing_state(rep_root):
    if not os.path.exists(rep_root):
        return False
    try:
        return any(True for _ in os.scandir(rep_root))
    except OSError:
        return True

