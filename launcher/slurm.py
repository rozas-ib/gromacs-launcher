import os
import re
import shutil
import subprocess


def shell_join(parts):
    return " ".join(str(part) for part in parts if str(part).strip())


def get_slurm_step_launcher(slurm_cfg, stage, step_name):
    stage_key = f"{step_name}_launcher_{stage}"
    if stage_key in slurm_cfg:
        return str(slurm_cfg.get(stage_key, "")).strip()
    shared_key = f"{step_name}_launcher"
    if shared_key in slurm_cfg:
        return str(slurm_cfg.get(shared_key, "")).strip()
    if step_name == "grompp":
        return "srun -n 1"
    if step_name == "mdrun":
        ntasks = slurm_cfg.get(f"mdrun_ntasks_{stage}")
        if ntasks:
            return f"srun -n {ntasks}"
        return "srun -n 1"
    return "srun -n 1"


def get_mdrun_extra_args(slurm_cfg, stage):
    stage_key = f"mdrun_args_{stage}"
    if stage_key in slurm_cfg:
        return str(slurm_cfg.get(stage_key, "")).strip()
    if "mdrun_args" in slurm_cfg:
        return str(slurm_cfg.get("mdrun_args", "")).strip()
    return ""


def build_sbatch_directives(slurm_cfg, stage, job_name, output_path, error_path):
    directives = [
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --account={slurm_cfg['account']}",
        f"#SBATCH --partition={slurm_cfg[f'partition_{stage}']}",
    ]

    qos = slurm_cfg.get(f"qos_{stage}")
    if qos:
        directives.append(f"#SBATCH --qos={qos}")

    mem = slurm_cfg.get(f"mem_{stage}")
    if mem:
        directives.append(f"#SBATCH --mem={mem}")

    directives.extend(
        [
            f"#SBATCH --time={slurm_cfg[f'time_{stage}']}",
            f"#SBATCH --nodes={slurm_cfg.get(f'nodes_{stage}', 1)}",
            f"#SBATCH --cpus-per-task={slurm_cfg.get(f'cpus_per_task_{stage}', 2)}",
            f"#SBATCH --output={output_path}",
            f"#SBATCH --error={error_path}",
        ]
    )

    ntasks = slurm_cfg.get(f"ntasks_{stage}")
    if ntasks:
        directives.insert(-3, f"#SBATCH --ntasks={ntasks}")
    else:
        directives.insert(-3, f"#SBATCH --ntasks-per-node={slurm_cfg.get(f'ntasks_per_node_{stage}', 56)}")
    return "\n".join(directives)


def write_setup_sh(path, cfg, abs_rep_path, aggregate_log_path):
    project_cfg, slurm_cfg = cfg["project_settings"], cfg["slurm_settings"]
    module_block = "\n".join(project_cfg.get("gmx_modules", []))
    analysis_script = os.path.join(os.getcwd(), project_cfg["scripts_dir"], project_cfg["density_analysis_script"])
    grompp_launcher = get_slurm_step_launcher(slurm_cfg, "setup", "grompp")
    mdrun_launcher = get_slurm_step_launcher(slurm_cfg, "setup", "mdrun")
    mdrun_extra_args = get_mdrun_extra_args(slurm_cfg, "setup")
    sbatch_block = build_sbatch_directives(
        slurm_cfg,
        "setup",
        f"S_{project_cfg['system_name']}",
        f"{abs_rep_path}/setup.out",
        f"{abs_rep_path}/setup.err",
    )

    script_content = f"""#!/bin/bash
{sbatch_block}

set -euo pipefail

{module_block}
source ~/miniconda3/bin/activate {project_cfg['conda_env']}

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

BASE="{abs_rep_path}"
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
if [ ! -f "min_out.gro" ]; then
    if ! {shell_join([grompp_launcher, "$GMX", "grompp -f min.mdp -c start.gro -p topol.top -o min.tpr"])} > $BASE/1_min/grompp_min.log 2>&1; then
        append_stage_log "Replica $(basename "$BASE") | 1_min/grompp_min.log" "$BASE/1_min/grompp_min.log"
        echo "ERROR: grompp failed for minimization. See $BASE/1_min/grompp_min.log"
        exit 1
    fi
    append_stage_log "Replica $(basename "$BASE") | 1_min/grompp_min.log" "$BASE/1_min/grompp_min.log"
    {shell_join([mdrun_launcher, "$GMX", "mdrun -deffnm min -c min_out.gro", mdrun_extra_args])}
else
    echo "Minimization already completed. Skipping."
fi

cd $BASE/2_press
if [ ! -f "press_out.gro" ]; then
    if ! {shell_join([grompp_launcher, "$GMX", "grompp -f press.mdp -c ../1_min/min_out.gro -p topol.top -o press.tpr"])} > $BASE/2_press/grompp_press.log 2>&1; then
        append_stage_log "Replica $(basename "$BASE") | 2_press/grompp_press.log" "$BASE/2_press/grompp_press.log"
        echo "ERROR: grompp failed for pressure equilibration. See $BASE/2_press/grompp_press.log"
        exit 1
    fi
    append_stage_log "Replica $(basename "$BASE") | 2_press/grompp_press.log" "$BASE/2_press/grompp_press.log"
    {shell_join([mdrun_launcher, "$GMX", "mdrun -deffnm press -c press_out.gro -rdd 1.2 -pin on", mdrun_extra_args])}
else
    echo "Pressure equilibration already completed. Skipping."
fi

cd $BASE/3_anneal
if [ ! -f "anneal_out.gro" ]; then
    if ! {shell_join([grompp_launcher, "$GMX", "grompp -f anneal.mdp -c ../2_press/press_out.gro -p topol.top -o anneal.tpr"])} > $BASE/3_anneal/grompp_anneal.log 2>&1; then
        append_stage_log "Replica $(basename "$BASE") | 3_anneal/grompp_anneal.log" "$BASE/3_anneal/grompp_anneal.log"
        echo "ERROR: grompp failed for annealing. See $BASE/3_anneal/grompp_anneal.log"
        exit 1
    fi
    append_stage_log "Replica $(basename "$BASE") | 3_anneal/grompp_anneal.log" "$BASE/3_anneal/grompp_anneal.log"
    {shell_join([mdrun_launcher, "$GMX", "mdrun -deffnm anneal -c anneal_out.gro -pin on", mdrun_extra_args])}
else
    echo "Annealing already completed. Skipping."
fi

cd $BASE/3_anneal
if [ ! -f "$BASE/4_prod/start.gro" ]; then
    python3 {analysis_script} anneal.edr --num_replicas 1 --time_crop {project_cfg['density_analysis_time_crop']}
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
    with open(os.path.join(path, "run_setup.sh"), "w", encoding="utf-8") as f:
        f.write(script_content)


def write_prod_sh(path, cfg, abs_rep_path, chunk_idx, aggregate_log_path):
    project_cfg, slurm_cfg = cfg["project_settings"], cfg["slurm_settings"]
    module_block = "\n".join(project_cfg.get("gmx_modules", []))
    prod_path = os.path.join(abs_rep_path, "4_prod")
    grompp_launcher = get_slurm_step_launcher(slurm_cfg, "prod", "grompp")
    mdrun_launcher = get_slurm_step_launcher(slurm_cfg, "prod", "mdrun")
    mdrun_extra_args = get_mdrun_extra_args(slurm_cfg, "prod")
    sbatch_block = build_sbatch_directives(
        slurm_cfg,
        "prod",
        f"P{chunk_idx}_{project_cfg['system_name']}",
        f"{prod_path}/chunk_{chunk_idx}.out",
        f"{prod_path}/chunk_{chunk_idx}.err",
    )
    script_content = f"""#!/bin/bash
{sbatch_block}

set -euo pipefail

{module_block}

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

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

cd {prod_path}
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
    {shell_join([mdrun_launcher, "$GMX", 'mdrun -v -deffnm "$DEFFNM" -s "$TPR_FILE" -cpi "$CPT_FILE" -append -maxh 71.7 -pin on', mdrun_extra_args])}
else
    if [ -n "$TPR_FILE" ]; then
        {shell_join([mdrun_launcher, "$GMX", 'mdrun -v -deffnm "$DEFFNM" -s "$TPR_FILE" -maxh 71.7 -pin on', mdrun_extra_args])}
    else
        if ! {shell_join([grompp_launcher, "$GMX", "grompp -f nvt_run.mdp -c start.gro -p topol.top -o nvt_run.tpr"])} > {prod_path}/grompp_prod.log 2>&1; then
            append_stage_log "Replica $(basename "$(dirname "$PWD")") | 4_prod/grompp_prod.log" "{prod_path}/grompp_prod.log"
            echo "ERROR: grompp failed for production. See {prod_path}/grompp_prod.log"
            exit 1
        fi
        append_stage_log "Replica $(basename "$(dirname "$PWD")") | 4_prod/grompp_prod.log" "{prod_path}/grompp_prod.log"
        {shell_join([mdrun_launcher, "$GMX", "mdrun -v -deffnm nvt_run -s nvt_run.tpr -maxh 71.7 -pin on", mdrun_extra_args])}
    fi
fi
"""
    with open(os.path.join(path, f"run_prod_{chunk_idx}.sh"), "w", encoding="utf-8") as f:
        f.write(script_content)


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
    setup_pattern = re.compile(rf"Replica {replica_idx}: submitted setup job with id ([^\s]+)")
    prod_pattern = re.compile(rf"Replica {replica_idx}: submitted production chunk \d+ with job id ([^\s]+)")
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
            ["squeue", "--noheader", "--format=%i|%T", "--jobs", ",".join(job_ids)],
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
            ["sacct", "-n", "-P", "-X", "-j", ",".join(job_ids), "--format=JobIDRaw,State"],
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


def submit_sbatch(script_path, dependency_job_id=None):
    cmd = ["sbatch", "--parsable"]
    if dependency_job_id:
        cmd.append(f"--dependency=afterany:{dependency_job_id}")
    cmd.append(script_path)
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except OSError as exc:
        raise RuntimeError(f"Failed to execute sbatch for '{script_path}': {exc}") from exc

    stdout = (res.stdout or "").strip()
    stderr = (res.stderr or "").strip()
    if res.returncode != 0:
        raise RuntimeError(
            f"sbatch failed for '{script_path}' with return code {res.returncode}. stderr: {stderr or '[empty]'}"
        )
    if not stdout:
        raise RuntimeError(
            f"sbatch returned no job id for '{script_path}'. stderr: {stderr or '[empty]'}"
        )
    return stdout
