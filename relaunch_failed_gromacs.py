#!/usr/bin/env python3

import argparse
import itertools
import os
import shutil
import sys

import launch_gromacs as lg


def classify_relaunch_candidate(row, rep_root):
    if not os.path.exists(rep_root):
        return "missing", False
    if row["slurm_status"] in {"queued", "running"}:
        return "active", False
    if row["prod_target_reached"] == "yes":
        return "completed", False
    return "failed_or_incomplete_idle", True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", nargs="?", default="config.toml")
    parser.add_argument(
        "--confirm",
        action="store_true",
        help="Actually relaunch the failed/incomplete idle replicas listed in the preview",
    )
    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"ERROR: config file '{args.config}' does not exist.")
        sys.exit(1)

    with open(args.config, "rb") as f:
        cfg = lg.tomllib.load(f)

    scr = cfg["screening"]
    sp = lg.normalize_species_config(cfg["species"])
    def to_list(val):
        return val if isinstance(val, list) else [val]
    group_defs, group_keys, active_species_keys = lg.discover_groups(cfg, sp)
    order = active_species_keys
    component_ratios = lg.normalize_component_ratios(scr, group_keys)
    force_field_sets = lg.normalize_force_field_sets(scr, order)
    master_combos = list(itertools.product(component_ratios, to_list(scr["target_temps"]), force_field_sets))
    output_root = lg.get_output_root(cfg)
    case_contexts = lg.build_case_contexts(master_combos, order, group_keys, output_root)
    progress_rows = lg.collect_progress_snapshot(case_contexts, cfg)

    rows_by_key = {
        (row["simulation"], int(row["replica"])): row
        for row in progress_rows
    }

    summary = {
        "missing": 0,
        "active": 0,
        "completed": 0,
        "failed_or_incomplete_idle": 0,
    }
    candidates = []
    for ctx in case_contexts:
        label = ctx["label"]
        for replica_idx in range(1, cfg["project_settings"]["num_replicas"] + 1):
            rep_root = os.path.join(ctx["sys_path"], f"rep_{replica_idx}")
            row = rows_by_key[(label, replica_idx)]
            classification, eligible = classify_relaunch_candidate(row, rep_root)
            summary[classification] += 1
            if eligible:
                candidates.append((ctx, replica_idx, row))

    print("Relaunch preview")
    print(f"Missing replicas: {summary['missing']}")
    print(f"Active replicas to skip: {summary['active']}")
    print(f"Completed replicas to skip: {summary['completed']}")
    print(f"Failed/incomplete idle replicas eligible for relaunch: {summary['failed_or_incomplete_idle']}")

    if candidates:
        print("\nCandidates:")
        for ctx, replica_idx, row in candidates:
            print(
                f"- {ctx['label']} R{replica_idx}: "
                f"slurm_status={row['slurm_status']}, "
                f"setup_start={row['analysis_start_gro']}, "
                f"prod_target_reached={row['prod_target_reached']}, "
                f"next_chunk={row['next_chunk_idx']}"
            )
    else:
        print("\nNo failed/incomplete idle replicas are eligible for relaunch.")

    if not args.confirm:
        print("\nNo submissions were made. Re-run with --confirm to relaunch the candidates above.")
        return

    if not candidates:
        print("\nNothing to relaunch.")
        return

    gmx_path = cfg["project_settings"].get("gmx_executable_path")
    if not gmx_path or not os.path.exists(gmx_path):
        gmx_path = shutil.which("gmx_mpi")
    if not gmx_path:
        print("ERROR: GMX not found")
        sys.exit(1)

    for ctx, replica_idx, row in candidates:
        ratio_entry = ctx["ratio_entry"]
        temp = ctx["temp"]
        ff_entry = ctx["ff_entry"]
        active_itps = ctx["active_itps"]
        label = ctx["label"]
        sys_path = ctx["sys_path"]
        rep_root = os.path.join(sys_path, f"rep_{replica_idx}")
        job_log, job_log_path = lg.make_job_logger(sys_path, label)
        grompp_log_path = lg.rebuild_replica_grompp_log(rep_root)
        job_log(f"Replica {replica_idx}: relaunch-failed mode selected this replica for recovery.")
        job_log(f"Replica {replica_idx}: aggregated grompp log path: {grompp_log_path}")

        counts, sizing_info = lg.resolve_system_size(
            cfg,
            sp,
            order,
            group_defs,
            group_keys,
            ratio_entry["weights"],
        )
        lg.prepare_replica_stage_inputs(rep_root, cfg, order, active_itps, sp, counts, temp, label)
        lg.build_initial_box_if_needed(
            rep_root,
            label,
            counts,
            order,
            sp,
            gmx_path,
            sizing_info["box_size_nm"],
            job_log,
            grompp_log_path,
            replica_idx,
        )

        prod_dir = os.path.join(rep_root, "4_prod")
        setup_done = lg.is_setup_complete(rep_root)
        has_start = os.path.exists(os.path.join(prod_dir, "start.gro"))
        has_checkpoint = os.path.exists(os.path.join(prod_dir, "nvt_run.cpt"))
        target_steps = lg.get_target_prod_steps(prod_dir, cfg)
        chunk_status = lg.parse_chunk_err_status(prod_dir, target_steps)

        prev = None
        if not setup_done:
            lg.write_setup_sh(rep_root, cfg, rep_root, grompp_log_path)
            try:
                sid = lg.submit_sbatch(os.path.join(rep_root, "run_setup.sh"))
            except RuntimeError as exc:
                job_log(f"Replica {replica_idx}: ERROR submitting setup job during relaunch: {exc}")
                print(f"ERROR submitting setup job for {label} R{replica_idx}: {exc}")
                sys.exit(1)
            print(f"Relaunched {label} R{replica_idx}: Setup {sid}")
            job_log(f"Replica {replica_idx}: submitted setup job with id {sid}")
            prev = sid
        else:
            job_log(f"Replica {replica_idx}: setup already complete during relaunch, setup submission skipped.")

        if has_start and not has_checkpoint:
            start_chunk_idx = lg.get_next_unused_chunk_idx(prod_dir)
            prod_jobs = cfg["project_settings"]["num_prod_chunks"]
        elif has_start and has_checkpoint:
            start_chunk_idx = max(
                chunk_status["next_chunk_idx"],
                lg.get_next_unused_chunk_idx(prod_dir),
            )
            prod_jobs = 1
        else:
            start_chunk_idx = lg.get_next_unused_chunk_idx(prod_dir)
            prod_jobs = cfg["project_settings"]["num_prod_chunks"]

        for chunk_idx in range(start_chunk_idx, start_chunk_idx + prod_jobs):
            lg.write_prod_sh(rep_root, cfg, rep_root, chunk_idx, grompp_log_path)
            try:
                prev = lg.submit_sbatch(
                    os.path.join(rep_root, f"run_prod_{chunk_idx}.sh"),
                    dependency_job_id=prev,
                )
            except RuntimeError as exc:
                job_log(f"Replica {replica_idx}: ERROR submitting production chunk {chunk_idx} during relaunch: {exc}")
                print(f"ERROR submitting production chunk {chunk_idx} for {label} R{replica_idx}: {exc}")
                sys.exit(1)
            job_log(f"Replica {replica_idx}: submitted production chunk {chunk_idx} with job id {prev}")

        print(
            f"Relaunched {label} R{replica_idx}: "
            f"{prod_jobs} production job(s) submitted (starting at chunk {start_chunk_idx})."
        )


if __name__ == "__main__":
    main()
