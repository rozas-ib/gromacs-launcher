import os
import sys

try:
    import tomllib
except ModuleNotFoundError:
    try:
        import tomli as tomllib
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "Python < 3.11 requires 'tomli' for TOML config parsing. Install it with: pip install tomli"
        ) from exc


DEFAULT_CONFIG = """# GROMACS launch configuration (TOML)
# Define every molecular species using [[species]] entries.
# The script derives `resname` and `atoms_per_molecule` automatically from each species GRO file.
# Define formulation groups under [groups.*] using species names and stoichiometric coefficients.
# System molar ratios are specified only between groups.
#
# Screened system compositions are defined in [[screening.component_ratios]]
# as weights over groups. Group internal stoichiometry stays in [groups.*].stoichiometric_ratio.
# Force-field choices (species -> ITP file mappings) are screened independently in
# [[screening.force_field]], allowing a Cartesian product of ratios x temperatures x FF sets.

[project_settings]
system_name = "LiFSI_DME_TOL"
output_root = "data"
num_replicas = 1
scripts_dir = "scripts"
density_analysis_script = "density_analysis_tool.py"
density_analysis_time_crop = "600:1000"
box_size_nm = 15.0
box_scale_factor = 1.5
num_prod_chunks = 4
gmx_modules = [
#    "module purge",
#    "module load bsc/1.0",
#    "module load intel/2024.2",
#    "module load impi/2021.13",
#    "module load mkl/2024.2",
    "module load oneapi/2024.2",
    "module load gromacs/2024.2"
]
conda_env = "mdanalysis"
global_ff_files = ["ffbonded.itp", "ffnonbonded.itp", "forcefield.itp"]
topology_forcefield_include = "forcefield.itp"
template_dir = "mdp_templates"
gmx_executable_path = "/apps/GPP/GROMACS/2024.2/INTEL24.2/IMPI21.13/bin/gmx_mpi"

[system_sizing]
# Options:
#   "target_atoms"      -> estimate counts from target atom count + ratios
#   "reference_component" -> keep a chosen top-level component count fixed (salt or solvent)

#mode = "target_atoms"
#target_atoms = 12000
#estimate_box_from_atoms = true # -> if set to true, insert-molecules will use this volume
#atom_number_density_atoms_per_nm3 = 85.0

mode = "reference_component"
# Keep one top-level component count fixed using the same names as in
# [groups.*] (example names are arbitrary):
reference_component_key = "LiFSI_salt"    # or "ZnTFSI2_salt", "Solvents", etc.
reference_component_count = 150
estimate_box_from_atoms = false # -> if set to true, insert-molecules will use this volume
atom_number_density_atoms_per_nm3 = 85.0
#
# If estimate_box_from_atoms = false, project_settings.box_size_nm is used.

[screening]
# Component ratio basis = arbitrary weights across groups.
# Example below corresponds to LiFSI_salt:DME_solv:TOL_solv = 1:2:2.
target_temps = [298.15]

[[screening.component_ratios]]
# Weights over groups, not necessarily normalized.
LiFSI_salt = 1.0
DME_solv = 2.0
TOL_solv = 2.0

# Add more entries to screen more total formulations (example with grouped solvents):
# [[screening.component_ratios]]
# LiFSI_salt = 1.0
# Solvents = 1.0

[[screening.force_field]]
name = "ff_set1"
cation = "li.itp"
anion = "fsi.itp"
solvent_1 = "dme_qforce_cm5_1.itp"
solvent_2 = "tol_qforce_cm5_1.itp"

[[screening.force_field]]
name = "ff_set2"
cation = "li.itp"
anion = "fsi.itp"
solvent_1 = "dme_qforce_cm5_09.itp"
solvent_2 = "tol_qforce_cm5_09.itp"

[slurm_settings]
account = "cice86"
partition_setup = "gpp"
partition_prod = "gpp"
qos_setup = "gp_debug"
qos_prod = "gp_resa"
time_setup = "02:00:00"
time_prod = "72:00:00"
nodes_setup = 1
ntasks_per_node_setup = 56
cpus_per_task_setup = 2
nodes_prod = 1
ntasks_per_node_prod = 56
cpus_per_task_prod = 2

[[species]]
# Species names are arbitrary. Group stoichiometric_ratio references these names.
name = "cation"
gro = "li.gro"

[[species]]
name = "anion"
gro = "fsi_freq.gro"

[[species]]
name = "solvent_1"
gro = "dme.gro"

[[species]]
name = "solvent_2"
gro = "tol.gro"

# Add more species as needed (example):
# [[species]]
# name = "acn"
# gro = "acn.gro"

[groups.LiFSI_salt]
# One unit of the LiFSI salt group contains:
stoichiometric_ratio = { cation = 1, anion = 1 }

[groups.DME_solv]
stoichiometric_ratio = { solvent_1 = 1 }

[groups.TOL_solv]
stoichiometric_ratio = { solvent_2 = 1 }

# Alternative grouped solvent blend example:
# [groups.Solvents]
# stoichiometric_ratio = { solvent_1 = 2, solvent_2 = 1 }

[simulation_settings.press]
ref_p = 10.0
ref_t = 10.0
nsteps = 150000
dt = 0.002

[simulation_settings.anneal]
ref_p = 1.0
temp_high = 500.0
annealing_time = "0 200 400 1000"
nsteps = 500000
dt = 0.002

[simulation_settings.prod]
nsteps = 500000
dt = 0.002
"""


def create_default_config(filename):
    """Create a template configuration TOML file if it doesn't exist."""
    with open(filename, "w", encoding="utf-8") as f:
        f.write(DEFAULT_CONFIG)
    print(f"\n[!] Template '{filename}' created.")
    sys.exit(0)


def load_config(config_path):
    if not os.path.exists(config_path):
        create_default_config(config_path)
    with open(config_path, "rb") as f:
        return tomllib.load(f)

