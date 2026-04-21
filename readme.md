## GROMACS MD Launcher

`launch_gromacs.py` prepares and submits screened GROMACS workflows on Slurm. It builds all combinations of:

- formulation ratios defined in `[[screening.component_ratios]]`
- target temperatures from `screening.target_temps`
- force-field sets from `[[screening.force_field]]`

For each case and replica, the script prepares stage folders, writes topology and tailored MDP files, builds the initial packed box with `gmx insert-molecules`, and submits setup and production jobs for initial submission only. Existing replica folders are left untouched by the initial launcher.

### Repository layout

The launcher expects this layout:

1. `launch_gromacs.py`
2. `inputs/`
3. `mdp_templates/`
4. `scripts/`

Required contents:

- `inputs/` must contain all species `.gro` files, all `.itp` files referenced by the selected force-field sets, and the global force-field files listed in `project_settings.global_ff_files`.
- `mdp_templates/` must contain `min.mdp`, `press.mdp`, `anneal.mdp`, and `nvt_run.mdp`.
- `scripts/` must contain the density analysis script named in `project_settings.density_analysis_script`.

### Configuration format

The launcher uses a TOML configuration file. If the file passed on the command line does not exist, the script creates a template TOML file and exits.

Default invocation:

```bash
python launch_gromacs.py
```

This creates `config.toml` if it is missing.

Custom config path:

```bash
python launch_gromacs.py my_config.toml
```

Recommended operational workflow:

```bash
python launch_gromacs.py config.toml
python launch_gromacs.py config.toml --track-progress
python relaunch_failed_gromacs.py config.toml --confirm
```

If you want a simulation campaign to be written somewhere other than `data/`, set:

```toml
[project_settings]
output_root = "data_ff_2_tailored"
```

This keeps different campaigns separated cleanly while preserving the same case-label logic inside the chosen root folder.

Dry-run mode prints the resolved case combinations and species counts without generating files, calling GROMACS, or submitting Slurm jobs:

```bash
python launch_gromacs.py config.toml --dry-run
```

Track-progress mode refreshes only the progress tables from the existing case folders, logs, and Slurm state. It does not submit jobs, rebuild input files, rewrite `grompp.log`, or interfere with running calculations:

```bash
python launch_gromacs.py config.toml --track-progress
```

Recovery mode is handled by a separate script. It previews the failed or incomplete idle replicas first, and only resubmits them if `--confirm` is present:

```bash
python relaunch_failed_gromacs.py config.toml
python relaunch_failed_gromacs.py config.toml --confirm
```

### Runtime environment

The launcher itself needs Python 3.11+ for `tomllib`, or Python < 3.11 with `tomli` installed.

Python packages required by this repository are listed in `requirements.txt`.

Recommended setup:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

This installs the Python dependencies used by:

- `launch_gromacs.py`
- `scripts/density_analysis_tool.py`

External tools are still required separately:

- GROMACS with `gmx_mpi` available either at `project_settings.gmx_executable_path` or on `PATH`
- Slurm commands such as `sbatch` and `srun`
- Any cluster-specific environment modules referenced in `project_settings.gmx_modules`
- Conda, if you use `project_settings.conda_env` inside the generated Slurm scripts

The submitted setup job script loads the runtime environment from the config:

- `project_settings.gmx_modules`: shell lines inserted near the top of the Slurm job script, typically `module load ...`
- `project_settings.conda_env`: Conda environment activated in the setup job
- `project_settings.gmx_executable_path`: absolute path to `gmx_mpi`; if it is missing or invalid during launcher execution, the script falls back to `gmx_mpi` from `PATH`

The generated Slurm scripts place all `#SBATCH` directives immediately after the shebang so scheduler options such as `--account`, partitions, and QoS are parsed correctly before the shell body starts.

There is no hardcoded `source ~/load_gmx24.sh` workflow in the launcher. Any module-loading commands should be declared in `project_settings.gmx_modules`.

### Configuration structure

The generated template shows the full schema. The main sections are:

- `[project_settings]`: workflow naming, replica count, paths, module commands, Conda environment, global force-field files, box defaults, and number of production chunks to submit
- `[project_settings].output_root`: top-level folder where all case folders and progress tables are written; defaults to `data`
- `[system_sizing]`: rules for turning formulation ratios into molecule counts
- `[screening]`: temperatures, component-ratio screening, and force-field screening
- `[slurm_settings]`: Slurm account, partitions, QoS, walltimes, and CPU layout
- `[[species]]`: species definitions using a `name` and a `.gro` file
- `[groups.*]`: top-level formulation groups with stoichiometric ratios over species
- `[simulation_settings.press]`, `[simulation_settings.anneal]`, `[simulation_settings.prod]`: stage-specific MDP replacements

### Species and groups

Species are declared with `[[species]]` entries:

```toml
[[species]]
name = "cation"
gro = "li.gro"
```

The launcher reads each GRO file and derives:

- the residue name used in `[ molecules ]`
- the atom count per molecule used by atom-based sizing

Groups are defined under `[groups.*]` and reference species names:

```toml
[groups.LiFSI_salt]
stoichiometric_ratio = { cation = 1, anion = 1 }
```

Formulation ratios are expressed between groups, not directly between species.

### Screening logic

Each launched case is the Cartesian product of:

1. one entry from `[[screening.component_ratios]]`
2. one value from `screening.target_temps`
3. one entry from `[[screening.force_field]]`

`[[screening.component_ratios]]` assigns weights to groups. Example:

```toml
[[screening.component_ratios]]
name = "ratio_1"
LiFSI_salt = 1.0
DME_solv = 2.0
TOL_solv = 2.0
```

`[[screening.force_field]]` maps every active species to an ITP file:

```toml
[[screening.force_field]]
name = "ff_set1"
cation = "li.itp"
anion = "fsi.itp"
solvent_1 = "dme_qforce_cm5_1.itp"
solvent_2 = "tol_qforce_cm5_1.itp"
```

Case directory names are generated from temperature, ratio-set name, force-field-set name, resolved group weights, and chosen ITP filenames.

### System sizing

The launcher supports two sizing modes in `[system_sizing]`:

1. `mode = "reference_component"`
2. `mode = "target_atoms"`

`reference_component` mode:

- fixes one top-level group count with `reference_component_key`
- uses `reference_component_count` as the anchor
- allocates the remaining group counts to preserve the requested ratio weights

`target_atoms` mode:

- estimates the total number of group units required to approach `target_atoms`
- uses the atom count derived from each species GRO file

Box size handling:

- if `estimate_box_from_atoms = false`, the script uses `project_settings.box_size_nm`
- if `estimate_box_from_atoms = true`, the script estimates the box volume from the total atom count and `atom_number_density_atoms_per_nm3`, then applies `project_settings.box_scale_factor`

### What the launcher writes

For every case and replica, the launcher creates:

- `<output_root>/<case_label>/launch.log`
- `<output_root>/<case_label>/rep_<N>/1_min`
- `<output_root>/<case_label>/rep_<N>/2_press`
- `<output_root>/<case_label>/rep_<N>/3_anneal`
- `<output_root>/<case_label>/rep_<N>/4_prod`

Inside each stage directory it copies:

- global force-field files
- the active species ITPs for the selected force-field set
- the species GRO files

It also writes:

- `topol.top` for every stage
- stage-specific MDP files copied or tailored from `mdp_templates/`
- `run_setup.sh`
- `run_prod_<chunk>.sh`
- `rep_<N>/grompp.log` in each replica directory as an aggregated log of `insert-molecules` and `grompp` stage logs for that replica
- `<output_root>/simulation_progress.md` and `<output_root>/simulation_progress.csv` as a stage-by-stage progress snapshot for all case/replica combinations in the current config

### Execution pipeline

For each replica, the launcher performs the following steps:

1. Create stage directories and copy/write all required input files.
2. Build `1_min/start.gro` locally with `gmx insert-molecules`.
3. Run pre-submission checks to detect existing setup outputs, production start files, checkpoints, and completed chunk logs.
4. Submit a setup Slurm job if setup is incomplete.
5. Submit one or more chained production jobs with `sbatch --dependency=afterany:...`.

If `sbatch` fails or returns an empty job id, the launcher now stops immediately and records the submission error in `launch.log` instead of silently continuing.

The setup job runs:

1. minimization in `1_min`
2. pressure equilibration in `2_press`
3. annealing in `3_anneal`
4. density analysis on `anneal.edr`, then moves the selected structure to `4_prod/start.gro`

The production job:

- reuses an existing `nvt_run.tpr` when present
- resumes from an existing checkpoint when the matching `.cpt` file exists
- otherwise runs `grompp` and starts a fresh production chunk
- exits immediately if `grompp` fails, leaving the previous logs intact for inspection
- safely handles the normal "no `.tpr` exists yet" case before the first production `grompp`

### Skip and resume behavior

The launcher is restart-aware and avoids redoing finished work when possible.

Initial submission mode:

- `launch_gromacs.py` is conservative and only submits replicas that do not already have an existing `rep_<N>` folder with contents
- if a replica folder already exists, the initial launcher skips it instead of trying to relaunch it
- failed or incomplete existing replicas should be handled with `relaunch_failed_gromacs.py`

Initial box build:

- if `rep_<N>/1_min/start.gro` already exists, the launcher skips rebuilding the minimization start structure

Setup stage:

- the setup job skips minimization if `1_min/min_out.gro` exists
- it skips pressure equilibration if `2_press/press_out.gro` exists
- it skips annealing if `3_anneal/anneal_out.gro` exists
- it skips density analysis if `4_prod/start.gro` already exists

Production stage:

- if chunk logs show the target production step count has already been reached, the replica is skipped entirely
- if setup is already complete, setup submission is skipped
- if `4_prod/start.gro` exists but no production checkpoint exists, the launcher submits the configured number of production chunks starting from the next unused chunk index
- if both `4_prod/start.gro` and `4_prod/nvt_run.cpt` exist, the launcher submits a single continuation chunk using a new chunk number so earlier `chunk_*.out/.err` files are never overwritten
- if `chunk_1.err` exists for an interrupted first production chunk and a checkpoint is present, a relaunch submits `chunk_2`, not `chunk_1`
- if `chunk_1` finished and `chunk_2` was interrupted, a relaunch submits `chunk_3`, preserving the logs from both previous chunks

### Logging and outputs

Key output files:

- `<output_root>/<case_label>/launch.log`: high-level launcher log for the case
- `<output_root>/<case_label>/rep_<N>/grompp.log`: aggregated replica-level log rebuilt from that replica's `insert-molecules.log` and `grompp_*.log` files, then appended by new setup/production `grompp` attempts for that replica
- `<output_root>/simulation_progress.md`: human-readable progress table with one row per simulation replica
- `<output_root>/simulation_progress.csv`: spreadsheet-friendly version of the same progress table
- `<output_root>/<case_label>/rep_<N>/1_min/insert-molecules.log`: packing log from `gmx insert-molecules`
- `<output_root>/<case_label>/rep_<N>/setup.out` and `setup.err`: setup Slurm stdout/stderr
- `<output_root>/<case_label>/rep_<N>/4_prod/chunk_<M>.out` and `chunk_<M>.err`: production Slurm stdout/stderr
- `<output_root>/<case_label>/rep_<N>/*/grompp_*.log`: `grompp` logs per stage

The chunk numbering is append-only for production resubmissions. If a production attempt stops early and you relaunch the workflow, the next submission uses a new chunk index instead of reusing the previous one. This preserves the full history of `chunk_*.out/.err` files.

`grompp.log` is intentionally different from `launch.log`:

- `launch.log` is append-only and keeps the full launcher decision history across reruns
- each replica has its own `rep_<N>/grompp.log`
- each replica `grompp.log` is overwritten at the start of each launcher execution
- after being overwritten, each replica `grompp.log` is repopulated from that replica's existing `insert-molecules.log` and `grompp_*.log` files
- new `grompp` attempts from the current setup or production submission are then appended to that replica-specific file

This means that if `insert-molecules` and minimization `grompp` succeeded, but pressure `grompp` failed for `rep_1`, a relaunch will regenerate `rep_1/grompp.log` with the stored `insert-molecules` and minimization `grompp` logs first, then append the new pressure `grompp` attempt from the resent setup job for that same replica.

The progress table is generated on every launcher run and summarizes each case/replica with columns such as:

- `slurm job id`
- `slurm status`
- `insert-molecules`
- `grompp min`
- `min run`
- `grompp press`
- `press run`
- `grompp anneal`
- `anneal run`
- `analysis/start.gro`
- `grompp prod`
- `prod cpt`
- `target steps reached`
- `max step seen`
- `target steps`
- `next chunk`

The launcher now also tries to populate the Slurm columns from the most recent known job id for each replica:

- it queries `squeue` first for live queue states such as `queued` or `running`
- if the job is no longer in `squeue`, it falls back to `sacct` when available to recover terminal states such as `completed`, `failed`, `cancelled`, or `timeout`
- if neither command is available on the current machine, the Slurm status falls back to `completed`, `unknown`, or `not_submitted` based on the files and recorded job ids

The rest of the table remains a filesystem snapshot. It reports what has already been created or completed according to the files and logs present in each replica directory.

If you want to refresh only this table while calculations are already running, use `--track-progress`. That mode is read-only with respect to the simulation workflow: it only observes the current files/logs and rewrites `<output_root>/simulation_progress.md` and `<output_root>/simulation_progress.csv`.

`relaunch_failed_gromacs.py` is the recovery entry point:

- it never creates new case folders for recovery
- it only considers already existing replica folders
- it skips replicas that are currently `queued` or `running`
- it skips replicas that already reached the production target steps
- without `--confirm`, it only prints a preview summary of the replicas eligible for recovery
- with `--confirm`, it resubmits only the failed or incomplete idle replicas

Common stage outputs:

- `1_min/min_out.gro`
- `2_press/press_out.gro`
- `3_anneal/anneal_out.gro`
- `4_prod/start.gro`
- `4_prod/nvt_run.tpr`
- `4_prod/nvt_run.cpt`

### Restart guide

The recommended way to resend an interrupted or failed calculation is:

```bash
python relaunch_failed_gromacs.py config.toml
python relaunch_failed_gromacs.py config.toml --confirm
```

The recovery script previews eligible replicas first, then only relaunches failed or incomplete idle replicas already present on disk when `--confirm` is provided.

Scenario 1: minimization, pressurization, and annealing finished, but production `grompp` failed

- setup is considered complete if `4_prod/start.gro` already exists, so setup is skipped
- production is resubmitted from the prepared production start structure
- the generated production script now exits immediately if `grompp` fails, so the failure is isolated to `grompp_prod.log` and does not continue into `mdrun`
- if previous chunk logs already exist, the retry uses the next unused chunk index instead of overwriting earlier `chunk_*.out/.err` files

Scenario 2: minimization and pressurization finished, but annealing `grompp` failed

- minimization and pressurization are skipped because their output GRO files already exist
- annealing is retried
- the generated setup script now exits immediately if annealing `grompp` fails, so the job stops cleanly and the error remains in `3_anneal/grompp_anneal.log`
- density analysis and production are only reached if annealing succeeds

Scenario 3: minimization, pressurization, and annealing finished, production started, but it did not reach the requested `nsteps`, and there is a checkpoint file

- on relaunch, the recovery script detects the checkpoint and treats the production as resumable
- the continuation submission uses a new chunk number instead of reusing the interrupted one
- for example, if `chunk_1.err` belongs to the interrupted attempt, the relaunch submits `chunk_2`
- the production script resumes with `mdrun -cpi ... -append`, so the trajectory/checkpoint continuation uses the existing production state while keeping a separate Slurm log file for the new attempt

Scenario 4: minimization, pressurization, and annealing finished, production finished one chunk, but a later chunk was interrupted

- on relaunch, the recovery script reads the existing `chunk_*.err` history to determine that production has not yet reached the target number of steps
- the continuation submission uses the next unused chunk number
- for example, if `chunk_1` finished and `chunk_2` was interrupted, the relaunch submits `chunk_3`
- this avoids overwriting `chunk_2.out` or `chunk_2.err` and keeps a complete record of the failed and resumed attempts

### Minimal usage workflow

1. Ensure `inputs/`, `mdp_templates/`, and `scripts/` contain the required files.
2. Generate a template config if needed:

```bash
python launch_gromacs.py config.toml
```

3. Edit `config.toml`.
4. Validate the case matrix without launching jobs:

```bash
python launch_gromacs.py config.toml --dry-run
```

5. Launch the workflow:

```bash
python launch_gromacs.py config.toml
```

6. Refresh the status tables without modifying the workflow:

```bash
python launch_gromacs.py config.toml --track-progress
```

7. Preview or confirm recovery of failed/incomplete idle replicas:

```bash
python relaunch_failed_gromacs.py config.toml
python relaunch_failed_gromacs.py config.toml --confirm
```

### Important parameters to review

Before launching, verify at least these fields:

- `project_settings.num_replicas`
- `project_settings.global_ff_files`
- `project_settings.template_dir`
- `project_settings.scripts_dir`
- `project_settings.density_analysis_script`
- `project_settings.density_analysis_time_crop`
- `project_settings.num_prod_chunks`
- `project_settings.gmx_modules`
- `project_settings.conda_env`
- `project_settings.gmx_executable_path`
- `system_sizing.mode`
- `system_sizing.reference_component_key`
- `system_sizing.reference_component_count`
- `system_sizing.target_atoms`
- `system_sizing.estimate_box_from_atoms`
- `screening.target_temps`
- `[[screening.component_ratios]]`
- `[[screening.force_field]]`
- `[[species]]`
- `[groups.*]`
- `[simulation_settings.press]`
- `[simulation_settings.anneal]`
- `[simulation_settings.prod]`
- `[slurm_settings]`
