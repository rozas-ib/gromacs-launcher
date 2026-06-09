import csv
import os
import re
from datetime import datetime

from .slurm import collect_slurm_status_map, parse_replica_job_ids


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


def is_setup_complete(rep_root):
    required = [
        os.path.join(rep_root, "1_min", "min_out.gro"),
        os.path.join(rep_root, "2_press", "press_out.gro"),
        os.path.join(rep_root, "3_anneal", "anneal_out.gro"),
        os.path.join(rep_root, "4_prod", "start.gro"),
    ]
    return all(os.path.exists(path) for path in required)


def read_mdp_value(mdp_path, key):
    if os.path.exists(mdp_path):
        try:
            with open(mdp_path, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    cleaned = line.split(";")[0].strip()
                    if "=" not in cleaned:
                        continue
                    item_key, val = [x.strip() for x in cleaned.split("=", 1)]
                    if item_key == key:
                        return val
        except OSError:
            pass
    return None


def read_mdp_int(mdp_path, key):
    value = read_mdp_value(mdp_path, key)
    if value is None:
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def read_mdp_float(mdp_path, key):
    value = read_mdp_value(mdp_path, key)
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def get_target_prod_steps(prod_dir, cfg):
    mdp_path = os.path.join(prod_dir, "nvt_run.mdp")
    nsteps = read_mdp_int(mdp_path, "nsteps")
    if nsteps is not None:
        return nsteps
    return int(cfg["simulation_settings"]["prod"]["nsteps"])


def parse_max_step_from_text(content):
    max_step = 0
    step_pattern = re.compile(r"\bstep\s+([0-9]+)\b", re.IGNORECASE)
    table_pattern = re.compile(r"^\s*([0-9]+)\s+[-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?(?:\s|$)")
    for match in step_pattern.finditer(content):
        max_step = max(max_step, int(match.group(1)))
    for line in content.splitlines():
        match = table_pattern.match(line)
        if match:
            max_step = max(max_step, int(match.group(1)))
    return max_step


def parse_max_step_from_files(paths):
    max_step = 0
    for path in paths:
        if not os.path.exists(path):
            continue
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                content = f.read()
        except (OSError, ValueError):
            continue
        max_step = max(max_step, parse_max_step_from_text(content))
    return max_step


def stage_time_progress(stage_dir, mdp_name, deffnm, output_gro, fallback_nsteps=None):
    mdp_path = os.path.join(stage_dir, mdp_name)
    target_steps = read_mdp_int(mdp_path, "nsteps")
    if target_steps is None and fallback_nsteps is not None:
        target_steps = int(fallback_nsteps)
    dt = read_mdp_float(mdp_path, "dt")
    completed = os.path.exists(os.path.join(stage_dir, output_gro))
    current_step = parse_max_step_from_files([
        os.path.join(stage_dir, f"{deffnm}.log"),
        os.path.join(stage_dir, f"{deffnm}.err"),
        os.path.join(stage_dir, f"{deffnm}.out"),
    ])
    if completed and target_steps is not None:
        current_step = max(current_step, target_steps)
    if target_steps is not None:
        current_step = min(current_step, target_steps)
    current_time_ps = "" if dt is None else current_step * dt
    target_time_ps = "" if dt is None or target_steps is None else target_steps * dt
    return {
        "current_step": current_step,
        "target_steps": "" if target_steps is None else target_steps,
        "current_time_ps": current_time_ps,
        "target_time_ps": target_time_ps,
    }


def parse_chunk_err_status(prod_dir, target_steps):
    completed_chunks = 0
    max_step = 0
    reached_target = False

    idx = 1
    last_seen_chunk_idx = 0
    while True:
        err_path = os.path.join(prod_dir, f"chunk_{idx}.err")
        if not os.path.exists(err_path):
            break
        last_seen_chunk_idx = idx
        try:
            with open(err_path, "r", encoding="utf-8", errors="ignore") as f:
                content = f.read()
        except OSError:
            break

        max_step = max(max_step, parse_max_step_from_text(content))
        if max_step >= target_steps:
            reached_target = True
            completed_chunks += 1
            break
        if ("Writing final coordinates." in content) or ("Finished mdrun" in content):
            completed_chunks += 1
            idx += 1
            continue
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


def production_time_progress(prod_dir, cfg, chunk_status):
    target_steps = get_target_prod_steps(prod_dir, cfg)
    dt = read_mdp_float(os.path.join(prod_dir, "nvt_run.mdp"), "dt")
    if dt is None:
        dt = float(cfg["simulation_settings"]["prod"]["dt"])
    log_step = parse_max_step_from_files([
        os.path.join(prod_dir, "nvt_run.log"),
        os.path.join(prod_dir, "nvt_run.err"),
        os.path.join(prod_dir, "nvt_run.out"),
    ])
    current_step = max(chunk_status["max_step"], log_step)
    if chunk_status["reached_target"]:
        current_step = max(current_step, target_steps)
    current_step = min(current_step, target_steps)
    return {
        "current_step": current_step,
        "target_steps": target_steps,
        "current_time_ps": current_step * dt,
        "target_time_ps": target_steps * dt,
    }


def status_flag(condition):
    return "yes" if condition else "no"


def build_progress_row(case_label, rep_idx, rep_root, cfg, launch_log_path, slurm_status_map):
    min_dir = os.path.join(rep_root, "1_min")
    press_dir = os.path.join(rep_root, "2_press")
    anneal_dir = os.path.join(rep_root, "3_anneal")
    prod_dir = os.path.join(rep_root, "4_prod")

    target_steps = get_target_prod_steps(prod_dir, cfg)
    chunk_status = parse_chunk_err_status(prod_dir, target_steps)
    min_progress = stage_time_progress(min_dir, "min.mdp", "min", "min_out.gro")
    press_progress = stage_time_progress(
        press_dir,
        "press.mdp",
        "press",
        "press_out.gro",
        fallback_nsteps=cfg["simulation_settings"]["press"]["nsteps"],
    )
    anneal_progress = stage_time_progress(
        anneal_dir,
        "anneal.mdp",
        "anneal",
        "anneal_out.gro",
        fallback_nsteps=cfg["simulation_settings"]["anneal"]["nsteps"],
    )
    prod_progress = production_time_progress(prod_dir, cfg, chunk_status)
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

    return {
        "simulation": case_label,
        "replica": rep_idx,
        "slurm_job_id": tracked_job_id or "",
        "slurm_status": slurm_status,
        "insert_molecules": status_flag(os.path.exists(os.path.join(min_dir, "start.gro"))),
        "grompp_min": status_flag(os.path.exists(os.path.join(min_dir, "min.tpr"))),
        "min_run": status_flag(os.path.exists(os.path.join(min_dir, "min_out.gro"))),
        "min_current_step": min_progress["current_step"],
        "min_target_steps": min_progress["target_steps"],
        "min_current_time_ps": min_progress["current_time_ps"],
        "min_target_time_ps": min_progress["target_time_ps"],
        "grompp_press": status_flag(os.path.exists(os.path.join(press_dir, "press.tpr"))),
        "press_run": status_flag(os.path.exists(os.path.join(press_dir, "press_out.gro"))),
        "press_current_step": press_progress["current_step"],
        "press_target_steps": press_progress["target_steps"],
        "press_current_time_ps": press_progress["current_time_ps"],
        "press_target_time_ps": press_progress["target_time_ps"],
        "grompp_anneal": status_flag(os.path.exists(os.path.join(anneal_dir, "anneal.tpr"))),
        "anneal_run": status_flag(os.path.exists(os.path.join(anneal_dir, "anneal_out.gro"))),
        "anneal_current_step": anneal_progress["current_step"],
        "anneal_target_steps": anneal_progress["target_steps"],
        "anneal_current_time_ps": anneal_progress["current_time_ps"],
        "anneal_target_time_ps": anneal_progress["target_time_ps"],
        "analysis_start_gro": status_flag(os.path.exists(os.path.join(prod_dir, "start.gro"))),
        "grompp_prod": status_flag(os.path.exists(os.path.join(prod_dir, "nvt_run.tpr"))),
        "prod_checkpoint": status_flag(os.path.exists(os.path.join(prod_dir, "nvt_run.cpt"))),
        "prod_target_reached": status_flag(chunk_status["reached_target"]),
        "prod_current_step": prod_progress["current_step"],
        "prod_target_steps": prod_progress["target_steps"],
        "prod_current_time_ps": prod_progress["current_time_ps"],
        "prod_target_time_ps": prod_progress["target_time_ps"],
        "max_step_seen": chunk_status["max_step"],
        "target_steps": target_steps,
        "next_chunk_idx": max(chunk_status["next_chunk_idx"], get_next_unused_chunk_idx(prod_dir)),
    }


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
        "min_current_step",
        "min_target_steps",
        "min_current_time_ps",
        "min_target_time_ps",
        "grompp_press",
        "press_run",
        "press_current_step",
        "press_target_steps",
        "press_current_time_ps",
        "press_target_time_ps",
        "grompp_anneal",
        "anneal_run",
        "anneal_current_step",
        "anneal_target_steps",
        "anneal_current_time_ps",
        "anneal_target_time_ps",
        "analysis_start_gro",
        "grompp_prod",
        "prod_checkpoint",
        "prod_target_reached",
        "prod_current_step",
        "prod_target_steps",
        "prod_current_time_ps",
        "prod_target_time_ps",
        "max_step_seen",
        "target_steps",
        "next_chunk_idx",
    ]

    with open(csv_path, "w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    header_labels = {
        "simulation": "simulation",
        "replica": "replica",
        "slurm_job_id": "slurm job id",
        "slurm_status": "slurm status",
        "insert_molecules": "insert-molecules",
        "grompp_min": "grompp min",
        "min_run": "min run",
        "min_current_step": "min step",
        "min_target_steps": "min target steps",
        "min_current_time_ps": "min time ps",
        "min_target_time_ps": "min target ps",
        "grompp_press": "grompp press",
        "press_run": "press run",
        "press_current_step": "press step",
        "press_target_steps": "press target steps",
        "press_current_time_ps": "press time ps",
        "press_target_time_ps": "press target ps",
        "grompp_anneal": "grompp anneal",
        "anneal_run": "anneal run",
        "anneal_current_step": "anneal step",
        "anneal_target_steps": "anneal target steps",
        "anneal_current_time_ps": "anneal time ps",
        "anneal_target_time_ps": "anneal target ps",
        "analysis_start_gro": "analysis/start.gro",
        "grompp_prod": "grompp prod",
        "prod_checkpoint": "prod cpt",
        "prod_target_reached": "target steps reached",
        "prod_current_step": "prod step",
        "prod_target_steps": "prod target steps",
        "prod_current_time_ps": "prod time ps",
        "prod_target_time_ps": "prod target ps",
        "max_step_seen": "max step seen",
        "target_steps": "target steps",
        "next_chunk_idx": "next chunk",
    }
    with open(md_path, "w", encoding="utf-8") as md_file:
        md_file.write("# Simulation Progress\n\n")
        md_file.write("Automatically generated by the launcher workflow.\n\n")
        md_file.write("| " + " | ".join(header_labels[h] for h in fieldnames) + " |\n")
        md_file.write("| " + " | ".join("---" for _ in fieldnames) + " |\n")
        for row in rows:
            md_file.write("| " + " | ".join(str(row[h]) for h in fieldnames) + " |\n")

    return csv_path, md_path


def collect_progress_snapshot(case_contexts, cfg):
    slurm_job_ids = []
    progress_row_refs = []
    for ctx in case_contexts:
        for replica_idx in range(1, cfg["project_settings"]["num_replicas"] + 1):
            rep_root = os.path.join(ctx.sys_path, f"rep_{replica_idx}")
            progress_row_refs.append((ctx.label, replica_idx, rep_root))
            job_ids = parse_replica_job_ids(ctx.job_log_path, replica_idx)
            if job_ids["setup_job_id"]:
                slurm_job_ids.append(job_ids["setup_job_id"])
            if job_ids["prod_job_id"]:
                slurm_job_ids.append(job_ids["prod_job_id"])

    slurm_status_map = collect_slurm_status_map(slurm_job_ids)
    log_path_by_label = {ctx.label: ctx.job_log_path for ctx in case_contexts}
    return [
        build_progress_row(label, replica_idx, rep_root, cfg, log_path_by_label[label], slurm_status_map)
        for label, replica_idx, rep_root in progress_row_refs
    ]
