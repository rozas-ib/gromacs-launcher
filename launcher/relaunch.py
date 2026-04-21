import argparse
import os
import shutil
import sys

from .case_matrix import build_case_contexts, build_master_combos, discover_groups, normalize_species_config, resolve_system_size
from .config import load_config
from .paths import get_output_root
from .slurm import submit_sbatch, write_prod_sh, write_setup_sh
from .status_reporting import collect_progress_snapshot, get_next_unused_chunk_idx, get_target_prod_steps, is_setup_complete, make_job_logger, parse_chunk_err_status
from .system_setup import build_initial_box_if_needed, prepare_replica_stage_inputs, rebuild_replica_grompp_log


def classify_relaunch_candidate(row, rep_root):
    if not os.path.exists(rep_root):
        return "missing_on_disk", False
    if row["slurm_status"] in {"queued", "running"}:
        return "currently_active_skip", False
    if row["prod_target_reached"] == "yes":
        return "already_completed_skip", False
    return "idle_incomplete_or_failed_relaunchable", True


def run_relaunch(config_path, confirm=False):
    if not os.path.exists(config_path):
        print(f"ERROR: config file '{config_path}' does not exist.")
        return 1

    cfg = load_config(config_path)
    species_cfg = normalize_species_config(cfg["species"])
    group_defs, group_keys, active_species_keys = discover_groups(cfg, species_cfg)
    order = active_species_keys
    master_combos = build_master_combos(cfg, group_keys, order)
    output_root = get_output_root(cfg)
    case_contexts = build_case_contexts(master_combos, order, group_keys, output_root)
    progress_rows = collect_progress_snapshot(case_contexts, cfg)

    rows_by_key = {(row["simulation"], int(row["replica"])): row for row in progress_rows}
    summary = {
        "missing_on_disk": 0,
        "currently_active_skip": 0,
        "already_completed_skip": 0,
        "idle_incomplete_or_failed_relaunchable": 0,
    }
    candidates = []

    for ctx in case_contexts:
        for replica_idx in range(1, cfg["project_settings"]["num_replicas"] + 1):
            rep_root = os.path.join(ctx.sys_path, f"rep_{replica_idx}")
            row = rows_by_key[(ctx.label, replica_idx)]
            classification, eligible = classify_relaunch_candidate(row, rep_root)
            summary[classification] += 1
            if eligible:
                candidates.append((ctx, replica_idx, row))

    print("Relaunch preview")
    print(f"Replica folders missing on disk: {summary['missing_on_disk']}")
    print(f"Replicas currently active in Slurm and skipped: {summary['currently_active_skip']}")
    print(f"Replicas already completed and skipped: {summary['already_completed_skip']}")
    print(f"Idle replicas that are incomplete or failed and can be relaunched: {summary['idle_incomplete_or_failed_relaunchable']}")

    if candidates:
        print("\nRelaunch candidates:")
        for ctx, replica_idx, row in candidates:
            print(
                f"- {ctx.label} R{replica_idx}: "
                f"slurm_status={row['slurm_status']}, "
                f"setup_ready={row['analysis_start_gro']}, "
                f"prod_target_reached={row['prod_target_reached']}, "
                f"next_chunk={row['next_chunk_idx']}"
            )
    else:
        print("\nNo idle failed/incomplete replicas are currently eligible for relaunch.")

    if not confirm:
        print("\nNo submissions were made. Re-run with --confirm to relaunch the candidates above.")
        return 0

    if not candidates:
        print("\nNothing to relaunch.")
        return 0

    gmx_path = cfg["project_settings"].get("gmx_executable_path")
    if not gmx_path or not os.path.exists(gmx_path):
        gmx_path = shutil.which("gmx_mpi")
    if not gmx_path:
        print("ERROR: GMX not found")
        return 1

    for ctx, replica_idx, _row in candidates:
        ratio_entry = ctx.ratio_entry
        temp = ctx.temp
        active_itps = ctx.active_itps
        label = ctx.label
        rep_root = os.path.join(ctx.sys_path, f"rep_{replica_idx}")
        job_log, _job_log_path = make_job_logger(ctx.sys_path, label)
        grompp_log_path = rebuild_replica_grompp_log(rep_root)
        job_log(f"Replica {replica_idx}: relaunch-failed mode selected this replica for recovery.")
        job_log(f"Replica {replica_idx}: aggregated grompp log path: {grompp_log_path}")

        counts, sizing_info = resolve_system_size(cfg, species_cfg, order, group_defs, group_keys, ratio_entry["weights"])
        prepare_replica_stage_inputs(rep_root, cfg, order, active_itps, species_cfg, counts, temp, label)
        build_initial_box_if_needed(
            rep_root,
            label,
            counts,
            order,
            species_cfg,
            gmx_path,
            sizing_info["box_size_nm"],
            job_log,
            grompp_log_path,
            replica_idx,
        )

        prod_dir = os.path.join(rep_root, "4_prod")
        setup_done = is_setup_complete(rep_root)
        has_start = os.path.exists(os.path.join(prod_dir, "start.gro"))
        has_checkpoint = os.path.exists(os.path.join(prod_dir, "nvt_run.cpt"))
        target_steps = get_target_prod_steps(prod_dir, cfg)
        chunk_status = parse_chunk_err_status(prod_dir, target_steps)

        prev = None
        if not setup_done:
            write_setup_sh(rep_root, cfg, rep_root, grompp_log_path)
            try:
                setup_job_id = submit_sbatch(os.path.join(rep_root, "run_setup.sh"))
            except RuntimeError as exc:
                job_log(f"Replica {replica_idx}: ERROR submitting setup job during relaunch: {exc}")
                print(f"ERROR submitting setup job for {label} R{replica_idx}: {exc}")
                return 1
            print(f"Relaunched {label} R{replica_idx}: Setup {setup_job_id}")
            job_log(f"Replica {replica_idx}: submitted setup job with id {setup_job_id}")
            prev = setup_job_id
        else:
            job_log(f"Replica {replica_idx}: setup already complete during relaunch, setup submission skipped.")

        if has_start and not has_checkpoint:
            start_chunk_idx = get_next_unused_chunk_idx(prod_dir)
            prod_jobs = cfg["project_settings"]["num_prod_chunks"]
        elif has_start and has_checkpoint:
            start_chunk_idx = max(chunk_status["next_chunk_idx"], get_next_unused_chunk_idx(prod_dir))
            prod_jobs = 1
        else:
            start_chunk_idx = get_next_unused_chunk_idx(prod_dir)
            prod_jobs = cfg["project_settings"]["num_prod_chunks"]

        for chunk_idx in range(start_chunk_idx, start_chunk_idx + prod_jobs):
            write_prod_sh(rep_root, cfg, rep_root, chunk_idx, grompp_log_path)
            try:
                prev = submit_sbatch(
                    os.path.join(rep_root, f"run_prod_{chunk_idx}.sh"),
                    dependency_job_id=prev,
                )
            except RuntimeError as exc:
                job_log(f"Replica {replica_idx}: ERROR submitting production chunk {chunk_idx} during relaunch: {exc}")
                print(f"ERROR submitting production chunk {chunk_idx} for {label} R{replica_idx}: {exc}")
                return 1
            job_log(f"Replica {replica_idx}: submitted production chunk {chunk_idx} with job id {prev}")

        print(
            f"Relaunched {label} R{replica_idx}: "
            f"{prod_jobs} production job(s) submitted (starting at chunk {start_chunk_idx})."
        )
    return 0


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", nargs="?", default="config.toml")
    parser.add_argument(
        "--confirm",
        action="store_true",
        help="Actually relaunch the failed/incomplete idle replicas listed in the preview",
    )
    args = parser.parse_args()
    sys.exit(run_relaunch(args.config, confirm=args.confirm))

