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


def get_target_prod_steps(prod_dir, cfg):
    mdp_path = os.path.join(prod_dir, "nvt_run.mdp")
    if os.path.exists(mdp_path):
        try:
            with open(mdp_path, "r", encoding="utf-8", errors="ignore") as f:
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
            with open(err_path, "r", encoding="utf-8", errors="ignore") as f:
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

    return {
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

