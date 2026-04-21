import argparse
import os
import shutil
import sys
from dataclasses import replace

from .case_matrix import (
    build_case_contexts,
    build_master_combos,
    discover_groups,
    format_final_molar_ratio,
    normalize_species_config,
    print_dry_run_case,
    resolve_system_size,
)
from .config import load_config
from .paths import get_output_root
from .slurm import collect_slurm_status_map, submit_sbatch, write_prod_sh, write_setup_sh
from .status_reporting import (
    build_progress_row,
    collect_progress_snapshot,
    get_next_unused_chunk_idx,
    get_target_prod_steps,
    is_setup_complete,
    make_job_logger,
    parse_chunk_err_status,
    write_progress_table,
)
from .system_setup import (
    build_initial_box_if_needed,
    prepare_replica_stage_inputs,
    rebuild_replica_grompp_log,
    replica_has_existing_state,
)
from .slurm import parse_replica_job_ids


def run_launch(config_path, dry_run=False, track_progress=False):
    if dry_run and track_progress:
        print("ERROR: --dry-run and --track-progress cannot be used together.")
        return 1

    cfg = load_config(config_path)
    species_cfg = normalize_species_config(cfg["species"])
    group_defs, group_keys, active_species_keys = discover_groups(cfg, species_cfg)
    order = active_species_keys
    master_combos = build_master_combos(cfg, group_keys, order)
    output_root = get_output_root(cfg)

    print(f"Setting {len(master_combos)} total variations...")
    case_contexts = build_case_contexts(master_combos, order, group_keys, output_root)

    if track_progress:
        progress_rows = collect_progress_snapshot(case_contexts, cfg)
        progress_csv_path, progress_md_path = write_progress_table(progress_rows, output_root)
        print(f"Track-progress mode: refreshed {progress_md_path} and {progress_csv_path}")
        return 0

    if not dry_run:
        gmx_path = cfg["project_settings"].get("gmx_executable_path")
        if not gmx_path or not os.path.exists(gmx_path):
            gmx_path = shutil.which("gmx_mpi")
        if not gmx_path:
            print("ERROR: GMX not found")
            return 1

    progress_row_refs = []
    slurm_job_ids = []
    updated_contexts = []

    for ctx in case_contexts:
        ratio_entry = ctx.ratio_entry
        temp = ctx.temp
        active_itps = ctx.active_itps
        ratio_name = ratio_entry["name"]
        ff_name = ctx.ff_entry["name"]
        label = ctx.label
        sys_path = ctx.sys_path
        job_log, job_log_path = make_job_logger(sys_path, label)
        job_log(f"Case details: ratio={ratio_name}, ff={ff_name}, temp={temp} K")
        job_log(f"Log file path: {job_log_path}")

        counts, sizing_info = resolve_system_size(cfg, species_cfg, order, group_defs, group_keys, ratio_entry["weights"])
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
        for replica_idx in range(1, cfg["project_settings"]["num_replicas"] + 1):
            job_ids = parse_replica_job_ids(job_log_path, replica_idx)
            if job_ids["setup_job_id"]:
                slurm_job_ids.append(job_ids["setup_job_id"])
            if job_ids["prod_job_id"]:
                slurm_job_ids.append(job_ids["prod_job_id"])

        updated_contexts.append(replace(ctx, job_log_path=job_log_path))

        if dry_run:
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

        for replica_idx in range(1, cfg["project_settings"]["num_replicas"] + 1):
            rep_root = os.path.join(sys_path, f"rep_{replica_idx}")
            if replica_has_existing_state(rep_root):
                msg = (
                    f"Replica {replica_idx}: existing replica folder detected at {rep_root}. "
                    "Initial-launch mode will not modify or resubmit existing replicas."
                )
                print(f"{label} R{replica_idx}: existing replica detected. Skipping initial submission.")
                job_log(msg)
                progress_row_refs.append((label, replica_idx, rep_root))
                continue

            grompp_log_path = rebuild_replica_grompp_log(rep_root)
            job_log("-" * 70)
            job_log(f"Replica {replica_idx}")
            job_log("-" * 70)
            job_log(f"Replica {replica_idx}: preparing stage folders and topology/MDP inputs in {rep_root}")
            job_log(f"Replica {replica_idx}: aggregated grompp log path: {grompp_log_path}")

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
            job_log("Pre-submission checks")
            job_log("-" * 70)
            job_log(
                f"Replica {replica_idx}: pre-submit checks -> setup_done={setup_done}, has_start={has_start}, "
                f"has_checkpoint={has_checkpoint}, target_steps={target_steps}, "
                f"max_step_seen={chunk_status['max_step']}, completed_chunks={chunk_status['completed_chunks']}, "
                f"next_chunk_idx={chunk_status['next_chunk_idx']}, "
                f"last_seen_chunk_idx={chunk_status['last_seen_chunk_idx']}"
            )

            if chunk_status["reached_target"]:
                msg = (
                    f"Skipping {label} R{replica_idx}: production reached target steps "
                    f"({chunk_status['max_step']} / {target_steps})."
                )
                print(msg)
                job_log(f"Replica {replica_idx}: {msg}")
                progress_row_refs.append((label, replica_idx, rep_root))
                continue

            prev = None
            if not setup_done:
                write_setup_sh(rep_root, cfg, rep_root, grompp_log_path)
                try:
                    setup_job_id = submit_sbatch(os.path.join(rep_root, "run_setup.sh"))
                except RuntimeError as exc:
                    job_log(f"Replica {replica_idx}: ERROR submitting setup job: {exc}")
                    print(f"ERROR submitting setup job for {label} R{replica_idx}: {exc}")
                    return 1
                slurm_job_ids.append(setup_job_id)
                print(f"Launched {label} R{replica_idx}: Setup {setup_job_id}")
                job_log(f"Replica {replica_idx}: submitted setup job with id {setup_job_id}")
                prev = setup_job_id
            else:
                print(f"{label} R{replica_idx}: setup already complete. Skipping setup submission.")
                job_log(f"Replica {replica_idx}: setup already complete, setup submission skipped.")

            if has_start and not has_checkpoint:
                start_chunk_idx = get_next_unused_chunk_idx(prod_dir)
                prod_jobs = cfg["project_settings"]["num_prod_chunks"]
            elif has_start and has_checkpoint:
                start_chunk_idx = max(chunk_status["next_chunk_idx"], get_next_unused_chunk_idx(prod_dir))
                prod_jobs = 1
            else:
                start_chunk_idx = get_next_unused_chunk_idx(prod_dir)
                prod_jobs = cfg["project_settings"]["num_prod_chunks"]
            job_log(
                f"Replica {replica_idx}: production policy -> start_chunk_idx={start_chunk_idx}, "
                f"prod_jobs_to_submit={prod_jobs}"
            )

            for chunk_idx in range(start_chunk_idx, start_chunk_idx + prod_jobs):
                write_prod_sh(rep_root, cfg, rep_root, chunk_idx, grompp_log_path)
                try:
                    prev = submit_sbatch(
                        os.path.join(rep_root, f"run_prod_{chunk_idx}.sh"),
                        dependency_job_id=prev,
                    )
                except RuntimeError as exc:
                    job_log(f"Replica {replica_idx}: ERROR submitting production chunk {chunk_idx}: {exc}")
                    print(f"ERROR submitting production chunk {chunk_idx} for {label} R{replica_idx}: {exc}")
                    return 1
                slurm_job_ids.append(prev)
                job_log(f"Replica {replica_idx}: submitted production chunk {chunk_idx} with job id {prev}")

            print(
                f"Launched {label} R{replica_idx}: {prod_jobs} production job(s) submitted "
                f"(starting at chunk {start_chunk_idx})."
            )
            job_log(
                f"Replica {replica_idx}: production submissions complete "
                f"(count={prod_jobs}, starting_chunk={start_chunk_idx})"
            )
            progress_row_refs.append((label, replica_idx, rep_root))

    if not progress_row_refs:
        for ctx in updated_contexts or case_contexts:
            for replica_idx in range(1, cfg["project_settings"]["num_replicas"] + 1):
                rep_root = os.path.join(ctx.sys_path, f"rep_{replica_idx}")
                progress_row_refs.append((ctx.label, replica_idx, rep_root))

    slurm_status_map = collect_slurm_status_map(slurm_job_ids)
    log_path_by_label = {ctx.label: ctx.job_log_path for ctx in (updated_contexts or case_contexts)}
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
    progress_csv_path, progress_md_path = write_progress_table(progress_rows, output_root)
    print(f"Simulation progress table written to {progress_md_path} and {progress_csv_path}")
    return 0


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", nargs="?", default="config.toml")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print computed species counts and exit without running GROMACS or submitting jobs",
    )
    parser.add_argument(
        "--track-progress",
        action="store_true",
        help="Do not submit or modify simulations; only refresh the progress tables from existing files/logs and Slurm state",
    )
    args = parser.parse_args()
    sys.exit(run_launch(args.config, dry_run=args.dry_run, track_progress=args.track_progress))
