###### GROMACS MD launcher Workflow - Documentation

This is a Python script to run GROMACS simulations workflows in slurm.
It automates stoichiometry calculations of electrolytes (ction, anion, solvent_1, solvent_2), initial box packing, mdp parameter setting and Slurm job chaining with automatic checkpoint resuming.

### 1. Directory Requirements

Ensure your root directory is organized as follows before execution:

1. `launch_gromacs.py`: The main script.
2. `scripts/`: Must contain the `density_analysis_tool.py`.
3. `mdp_templates/`: Must contain `min.mdp`, `press.mdp`, `anneal.mdp`, and `nvt_run.mdp`.
4. `inputs/`: Must contain all `.itp` and `.gro` files for your species, plus global forcefield files (e.g., `ffbonded.itp`, `ffnonbonded.itp`, `forcefield.itp`).

### 2. Getting Started - Mare Nostrum 5 specific

1. **Prepare the Environment**: Load GROMACS and activate your Conda environment.

   ```bash
   source ~/load_gmx24.sh
   conda activate `env-name` 
   ```

2. **Generate a Template Config**: Run the script once to create a default JSON file.

   ```bash
   python launch_gromacs.py config.json
   ```

3. **Configure your system**: Edit `config.json` to define your concentrations, temperature range, model list and mdp file settings.

4. **Launch Workflow**: Execute the script with your configuration file.

   ```bash
   python launch_gromacs.py config.json
   ```

### 3. Core Logic & Bundling

The script uses a specific **"Bundle-then-Cross"** logic to ensure charge models (introduced through a list in itp files) are handled correctly without incorrect mixing:

1. **Salt Bundle (Lock-Step)**: The script pairs the "cation" ITP list and "anion" ITP list by index. (e.g., Cation Index 0 always stays with Anion Index 0).

2. **Solvent Bundle (Lock-Step)**: The script pairs "solvent_1" and "solvent_2" ITP lists by index.

3. **Cartesian Product**: The script calculates every possible combination between the **Salt Bundle** and the **Solvent Bundle**, crossed with every **Concentration Ratio** and **Temperature**.

4. **Strict Order**: To prevent GROMACS topology errors, the script enforces a strict ordering for both box insertion and topology writing: `Cation -> Anion -> Solvent 1 -> Solvent 2`.


### 4. The Simulation Pipeline

For every unique variation and replica, the script manages the following stages:

1. **Box Building**: Performed locally on the login node using `gmx insert-molecules`. A unique random seed is used for every replica. Logs are saved to `insert-molecules.log`.

2. **Setup Job (S_...)**: A Slurm job that executes Minimization, Pressurization (NPT), and Annealing (NPT, Temperature Ramp) sequentially.

3. **Density Analysis**: Integrated at the end of the Setup job. It automatically finds the frame closest to the average density of the equilibrated system and moves it to the production folder.

4. **Production Chunks (P1, P2...)**: Chained Slurm jobs using `--dependency=afterany`. Each chunk detects existing `.cpt` files to resume the simulation seamlessly until the total step count is reached.

### 5. Output Organization & Logging

All simulation data is stored in the `data/` directory. Logs are redirected to subfolders to keep the root directory clean:

1. **Box Building**: `data/[Variation]/rep_N/1_min/insert-molecules.log`

2. **Slurm Setup Log**: `data/[Variation]/rep_N/setup.out`

3. **Slurm Production Logs**: `data/[Variation]/rep_N/4_prod/chunk_N.out`

4. **GROMACS grompp Logs**: Found inside each stage folder (e.g., `grompp_min.log`, `grompp_anneal.log`).

5. **Density Plot**: `data/[Variation]/rep_N/3_anneal/density_analysis.png`.

### 6. Critical Configuration Parameters

-`num_replicas`: Number of replicas per system configuration to run.

-`density_analysis_time_crop`: Time range, in ps, the density_analysis_scripr will use to calculate the average density.

-`n_solvent_total`: The total count of Solvent 1 + Solvent 2 molecules. Salt count is calculated relative to this.

-`box_size_nm`: The initial box length. Ensure this is large enough to prevent packing failures.

-`num_prod_chunks`: The number of consecutive 72-hour jobs to submit.

-`gmx_executable_path`: The absolute path to the `gmx_mpi` binary.


