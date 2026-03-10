## GROMACS MD Launcher

`launch_gromacs.py` prepares and submits screened GROMACS workflows on Slurm. It builds all combinations of:

- formulation ratios defined in `[[screening.component_ratios]]`
- target temperatures from `screening.target_temps`
- force-field sets from `[[screening.force_field]]`

For each case and replica, the script prepares stage folders, writes topology and tailored MDP files, builds the initial packed box with `gmx insert-molecules`, submits setup and production jobs, and resumes existing runs when possible.

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

Dry-run mode prints the resolved case combinations and species counts without generating files, calling GROMACS, or submitting Slurm jobs:

```bash
python launch_gromacs.py config.toml --dry-run
```

### Runtime environment

The launcher itself needs Python 3.11+ for `tomllib`, or Python < 3.11 with `tomli` installed.

The submitted setup job script loads the runtime environment from the config:

- `project_settings.gmx_modules`: shell lines inserted near the top of the Slurm job script, typically `module load ...`
- `project_settings.conda_env`: Conda environment activated in the setup job
- `project_settings.gmx_executable_path`: absolute path to `gmx_mpi`; if it is missing or invalid during launcher execution, the script falls back to `gmx_mpi` from `PATH`

There is no hardcoded `source ~/load_gmx24.sh` workflow in the launcher. Any module-loading commands should be declared in `project_settings.gmx_modules`.

### Configuration structure

The generated template shows the full schema. The main sections are:

- `[project_settings]`: workflow naming, replica count, paths, module commands, Conda environment, global force-field files, box defaults, and number of production chunks to submit
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

- `data/<case_label>/launch.log`
- `data/<case_label>/rep_<N>/1_min`
- `data/<case_label>/rep_<N>/2_press`
- `data/<case_label>/rep_<N>/3_anneal`
- `data/<case_label>/rep_<N>/4_prod`

Inside each stage directory it copies:

- global force-field files
- the active species ITPs for the selected force-field set
- the species GRO files

It also writes:

- `topol.top` for every stage
- stage-specific MDP files copied or tailored from `mdp_templates/`
- `run_setup.sh`
- `run_prod_<chunk>.sh`

### Execution pipeline

For each replica, the launcher performs the following steps:

1. Create stage directories and copy/write all required input files.
2. Build `1_min/start.gro` locally with `gmx insert-molecules`.
3. Run pre-submission checks to detect existing setup outputs, production start files, checkpoints, and completed chunk logs.
4. Submit a setup Slurm job if setup is incomplete.
5. Submit one or more chained production jobs with `sbatch --dependency=afterany:...`.

The setup job runs:

1. minimization in `1_min`
2. pressure equilibration in `2_press`
3. annealing in `3_anneal`
4. density analysis on `anneal.edr`, then moves the selected structure to `4_prod/start.gro`

The production job:

- reuses an existing `nvt_run.tpr` when present
- resumes from an existing checkpoint when the matching `.cpt` file exists
- otherwise runs `grompp` and starts a fresh production chunk

### Skip and resume behavior

The launcher is restart-aware and avoids redoing finished work when possible.

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
- if `4_prod/start.gro` exists but no production checkpoint exists, the launcher submits the configured number of production chunks starting from chunk 1
- if both `4_prod/start.gro` and `4_prod/nvt_run.cpt` exist, the launcher submits a single continuation chunk starting from the next missing chunk index inferred from existing `chunk_*.err` files

### Logging and outputs

Key output files:

- `data/<case_label>/launch.log`: high-level launcher log for the case
- `data/<case_label>/rep_<N>/1_min/insert-molecules.log`: packing log from `gmx insert-molecules`
- `data/<case_label>/rep_<N>/setup.out` and `setup.err`: setup Slurm stdout/stderr
- `data/<case_label>/rep_<N>/4_prod/chunk_<M>.out` and `chunk_<M>.err`: production Slurm stdout/stderr
- `data/<case_label>/rep_<N>/*/grompp_*.log`: `grompp` logs per stage

Common stage outputs:

- `1_min/min_out.gro`
- `2_press/press_out.gro`
- `3_anneal/anneal_out.gro`
- `4_prod/start.gro`
- `4_prod/nvt_run.tpr`
- `4_prod/nvt_run.cpt`

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
