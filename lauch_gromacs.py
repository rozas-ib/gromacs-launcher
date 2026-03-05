import os, shutil, subprocess, itertools, sys, argparse

try:
    import tomllib
except ModuleNotFoundError:
    try:
        import tomli as tomllib
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "Python < 3.11 requires 'tomli' for TOML config parsing. Install it with: pip install tomli"
        ) from exc

#This script 
#
#


def create_default_config(filename):
    """Creates a template configuration TOML file if it doesn't exist."""
    default_cfg = """# GROMACS launch configuration (TOML)
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
num_replicas = 1
scripts_dir = "scripts"
density_analysis_script = "density_analysis_tool.py"
density_analysis_time_crop = "600:1000"
box_size_nm = 15.0
box_scale_factor = 1.5
num_prod_chunks = 4
gmx_modules = [
    "module purge",
    "module load bsc/1.0",
    "module load intel/2024.2",
    "module load impi/2021.13",
    "module load mkl/2024.2",
    "module load gromacs/2024.2"
]
conda_env = "mdanalysis"
global_ff_files = ["ffbonded.itp", "ffnonbonded.itp", "forcefield.itp"]
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
ntasks_setup = 24
ntasks_prod = 48

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
    with open(filename, "w") as f:
        f.write(default_cfg)
    print(f"\n[!] Template '{filename}' created.")
    sys.exit(0)

def modify_mdp(template_path, output_path, replacements):
    """Replaces MDP keys with tailored values."""
    with open(template_path, 'r') as f: lines = f.readlines()
    new_lines = []
    for line in lines:
        cleaned = line.split(';')[0].strip()
        if "=" in cleaned:
            key = cleaned.split('=')[0].strip()
            if key in replacements:
                new_lines.append(f"{key:<24}= {replacements[key]} ; Tailored\n")
                continue
        new_lines.append(line)
    with open(output_path, 'w') as f: f.writelines(new_lines)

def read_gro_species_metadata(gro_path):
    with open(gro_path, "r") as f:
        lines = f.readlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid GRO file '{gro_path}': too few lines")
    try:
        atom_count = int(lines[1].strip())
    except ValueError as exc:
        raise ValueError(f"Invalid GRO file '{gro_path}': line 2 must be atom count") from exc
    atom_lines = lines[2:2 + atom_count]
    if len(atom_lines) != atom_count:
        raise ValueError(f"Invalid GRO file '{gro_path}': expected {atom_count} atom lines")
    residue_names = {line[5:10].strip() for line in atom_lines if len(line) >= 10}
    if not residue_names:
        raise ValueError(f"Invalid GRO file '{gro_path}': could not parse residue name")
    if len(residue_names) != 1:
        raise ValueError(
            f"GRO file '{gro_path}' contains multiple residue names ({', '.join(sorted(residue_names))}); "
            "cannot auto-derive a single species name for topology output"
        )
    return {"resname": next(iter(residue_names)), "atoms_per_molecule": atom_count}

def normalize_species_config(raw_species, inputs_dir="inputs"):
    if not isinstance(raw_species, list) or not raw_species:
        raise ValueError(
            "species must be defined as a non-empty array of tables using [[species]] entries"
        )
    species_cfg = {}
    for entry in raw_species:
        if not isinstance(entry, dict):
            raise ValueError("Each [[species]] entry must be a table")
        name = entry.get("name")
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Each [[species]] entry must define a non-empty string 'name'")
        key = name.strip()
        if key in species_cfg:
            raise ValueError(f"Duplicate species name '{key}' in [[species]] entries")
        normalized = dict(entry)
        normalized["name"] = key
        gro_name = normalized.get("gro")
        if not isinstance(gro_name, str) or not gro_name.strip():
            raise ValueError(f"Species '{key}' must define a non-empty 'gro' filename")
        gro_name = gro_name.strip()
        normalized["gro"] = gro_name
        gro_meta = read_gro_species_metadata(os.path.join(inputs_dir, gro_name))
        normalized["resname"] = gro_meta["resname"]
        normalized["atoms_per_molecule"] = gro_meta["atoms_per_molecule"]
        species_cfg[key] = normalized
    return species_cfg

def discover_groups(cfg, species_cfg):
    groups_cfg = cfg.get("groups")
    if not groups_cfg or not isinstance(groups_cfg, dict):
        raise ValueError("[groups] must be defined as a table of named groups")
    group_defs = {}
    active_species_keys = []
    for group_key, group_entry in groups_cfg.items():
        if not isinstance(group_entry, dict):
            raise ValueError(f"groups.{group_key} must be a table")
        stoich = group_entry.get("stoichiometric_ratio")
        if not isinstance(stoich, dict) or not stoich:
            raise ValueError(f"groups.{group_key}.stoichiometric_ratio must be a non-empty table")
        norm = {}
        for sp_key, coeff in stoich.items():
            if sp_key not in species_cfg:
                raise KeyError(f"groups.{group_key}.stoichiometric_ratio references unknown species '{sp_key}'")
            coeff_i = int(coeff)
            if coeff_i <= 0:
                raise ValueError(f"groups.{group_key}.stoichiometric_ratio.{sp_key} must be > 0")
            norm[sp_key] = coeff_i
            if sp_key not in active_species_keys:
                active_species_keys.append(sp_key)
        group_defs[group_key] = {"stoichiometric_ratio": norm}
    if not active_species_keys:
        raise ValueError("[groups] must reference at least one species")
    return group_defs, list(group_defs.keys()), active_species_keys

def normalize_component_ratios(scr, component_keys):
    comps = scr.get("component_ratios")
    if not isinstance(comps, list) or not comps:
        raise KeyError(
            "Missing screening.component_ratios. "
            "Use [[screening.component_ratios]] with weights for groups."
        )
    out = []
    for idx, comp in enumerate(comps, start=1):
        if not isinstance(comp, dict):
            raise ValueError("Each screening.component_ratios entry must be a table")
        entry_name = comp.get("name", f"ratio_{idx}")
        unknown = [k for k in comp if k not in component_keys and k != "name"]
        if unknown:
            raise KeyError(f"Unknown component keys in ratio: {', '.join(unknown)}")
        weights = {k: float(comp.get(k, 0.0)) for k in component_keys}
        if sum(weights.values()) <= 0:
            raise ValueError("Each component ratio entry must have positive total weight")
        out.append({"name": str(entry_name), "weights": weights})
    return out

def normalize_force_field_sets(scr, species_keys):
    ff_sets = scr.get("force_field")
    if not isinstance(ff_sets, list) or not ff_sets:
        raise KeyError(
            "Missing screening.force_field. "
            "Use [[screening.force_field]] entries mapping species names to ITP files."
        )
    out = []
    for idx, entry in enumerate(ff_sets, start=1):
        if not isinstance(entry, dict):
            raise ValueError("Each screening.force_field entry must be a table")
        name = str(entry.get("name", f"ff_{idx}"))
        unknown = [k for k in entry if k not in species_keys and k != "name"]
        if unknown:
            raise KeyError(f"Unknown species keys in screening.force_field '{name}': {', '.join(unknown)}")
        missing = [k for k in species_keys if k not in entry]
        if missing:
            raise KeyError(f"Missing species ITP mappings in screening.force_field '{name}': {', '.join(missing)}")
        species_itps = {}
        for k in species_keys:
            v = entry[k]
            if not isinstance(v, str) or not v.strip():
                raise ValueError(f"screening.force_field '{name}' entry for species '{k}' must be a non-empty string")
            species_itps[k] = v.strip()
        out.append({"name": name, "species_itps": species_itps})
    return out

def allocate_weighted_counts(n_total, composition):
    total_weight = sum(float(v) for v in composition.values())
    if total_weight <= 0:
        raise ValueError("Composition total weight must be positive")
    raw = {k: (n_total * float(v) / total_weight) for k, v in composition.items()}
    counts = {k: int(v) for k, v in raw.items()}
    remainder = n_total - sum(counts.values())
    if remainder > 0:
        rank = sorted(raw, key=lambda k: (raw[k] - counts[k]), reverse=True)
        for k in rank[:remainder]:
            counts[k] += 1
    return counts

def accumulate_species_counts_from_groups(group_counts, group_defs, active_species_keys):
    species_counts = {k: 0 for k in active_species_keys}
    for group_key, n_units in group_counts.items():
        for sp_key, coeff in group_defs[group_key]["stoichiometric_ratio"].items():
            species_counts[sp_key] += int(n_units) * int(coeff)
    return species_counts

def get_atoms_per_molecule(species_cfg, species_key):
    atoms = species_cfg[species_key].get("atoms_per_molecule")
    if atoms is None:
        raise KeyError(
            f"Missing species.{species_key}.atoms_per_molecule required for system_sizing.mode='target_atoms'"
        )
    atoms = int(atoms)
    if atoms <= 0:
        raise ValueError(f"species.{species_key}.atoms_per_molecule must be > 0")
    return atoms

def estimate_box_size_nm(total_atoms, project_cfg, sizing_cfg):
    if sizing_cfg.get("estimate_box_from_atoms", False):
        atom_density = float(sizing_cfg.get("atom_number_density_atoms_per_nm3", 0.0))
        if atom_density <= 0:
            raise ValueError("system_sizing.atom_number_density_atoms_per_nm3 must be > 0 when estimate_box_from_atoms = true")
        volume_nm3 = float(total_atoms) / atom_density
        scaled_volume = project_cfg.get("box_scale_factor") * volume_nm3
        return (float(project_cfg.get("box_scale_factor",1.0)) * volume_nm3) ** (1.0 / 3.0), "estimated_and_scaled_from_atom_density"
    return float(project_cfg["box_size_nm"]), "project_settings.box_size_nm"

def resolve_system_size(cfg, species_cfg, order, group_defs, group_keys, component_ratio):
    project_cfg = cfg["project_settings"]
    sizing_cfg = cfg.get("system_sizing", {})
    mode = sizing_cfg.get("mode", "target_atoms")
    total_ratio_weight = sum(float(v) for v in component_ratio.values())
    if total_ratio_weight <= 0:
        raise ValueError("Component ratio must contain positive total weight")

    if mode == "reference_component":
        ref_component_key = sizing_cfg["reference_component_key"]
        ref_component_count = int(sizing_cfg["reference_component_count"])
        if ref_component_key not in group_keys:
            raise KeyError(
                f"system_sizing.reference_component_key '{ref_component_key}' is not a valid component. "
                f"Use one of: {', '.join(group_keys)}"
            )
        if ref_component_count < 0:
            raise ValueError("system_sizing.reference_component_count must be >= 0")
        ref_fraction = float(component_ratio[ref_component_key]) / total_ratio_weight
        if ref_fraction <= 0:
            raise ValueError(
                f"Component '{ref_component_key}' has zero weight in this screening.component_ratios entry"
            )
        n_total_components = int(round(ref_component_count / ref_fraction))
        component_counts = allocate_weighted_counts(n_total_components, component_ratio)
    elif mode == "target_atoms":
        target_atoms = int(sizing_cfg["target_atoms"])
        if target_atoms <= 0:
            raise ValueError("system_sizing.target_atoms must be > 0")
        avg_component_atoms = 0.0
        for group_key in group_keys:
            frac = float(component_ratio[group_key]) / total_ratio_weight
            group_atoms = sum(
                int(coeff) * get_atoms_per_molecule(species_cfg, sp_key)
                for sp_key, coeff in group_defs[group_key]["stoichiometric_ratio"].items()
            )
            avg_component_atoms += frac * group_atoms
        if avg_component_atoms <= 0:
            raise ValueError("Invalid sizing denominator. Check component ratios and species atom counts.")
        n_total_components = max(1, int(round(target_atoms / avg_component_atoms)))
        component_counts = allocate_weighted_counts(n_total_components, component_ratio)
    else:
        raise ValueError("system_sizing.mode must be one of: target_atoms, reference_component")
    group_counts = {k: int(component_counts[k]) for k in group_keys}
    counts = accumulate_species_counts_from_groups(group_counts, group_defs, order)

    total_atoms = None
    if mode == "target_atoms" or sizing_cfg.get("estimate_box_from_atoms", False):
        total_atoms = sum(counts[k] * get_atoms_per_molecule(species_cfg, k) for k in order)

    box_size_nm, box_source = (
        estimate_box_size_nm(total_atoms, project_cfg, sizing_cfg)
        if total_atoms is not None else
        (float(project_cfg["box_size_nm"]), "project_settings.box_size_nm")
    )

    return counts, {
        "mode": mode,
        "n_total_components": sum(group_counts.values()),
        "total_atoms": total_atoms,
        "box_size_nm": box_size_nm,
        "box_source": box_source,
        "component_counts": group_counts,
    }

def fmt_num(value):
    value = float(value)
    return str(int(value)) if value.is_integer() else f"{value:.4g}"

def build_case_label(temp, ratio_name, ff_name, component_ratio, active_itps, species_order, group_keys):
    comp_mix_label = "_".join(f"{k}{fmt_num(component_ratio[k])}" for k in group_keys)
    model_keys = species_order
    itp_label = "_".join(active_itps[k].split('.')[0] for k in model_keys)
    return f"T{fmt_num(temp)}_{ratio_name}_{ff_name}_{comp_mix_label}_{itp_label}"

def print_dry_run_case(label, temp, ratio_name, ff_name, counts, order, group_keys, component_ratio, active_itps, sizing_info):
    total_system = sum(counts[k] for k in order)
    component_weights = ", ".join(f"{k}={fmt_num(component_ratio[k])}" for k in group_keys)
    species_counts = ", ".join(f"{k}={counts[k]}" for k in order)
    print(f"{'='*100}\n[DRY-RUN] {label}\n{'='*100}")
    print(f"Temp={fmt_num(temp)} K\n Ratio={ratio_name}\n FF={ff_name}\n Component_ratio=({component_weights})")
    details = (
        f"Sizing_mode={sizing_info['mode']} -->  Total_groups={sizing_info['n_total_components']}, "
        f"box_nm={sizing_info['box_size_nm']:.3f} ({sizing_info['box_source']})"
    )
    print(details)
    print("  itp_map: " + ", ".join(f"{k}={active_itps[k]}" for k in order))
    component_counts = sizing_info.get("component_counts")
    if component_counts:
        print("  groups: " + ", ".join(f"{k}={component_counts[k]}" for k in group_keys))
    print(f"  counts: {species_counts}")
    totals_line = f"  totals: system={total_system}"
    if sizing_info["total_atoms"] is not None:
        totals_line += f", atoms={sizing_info['total_atoms']}"
    print(totals_line)

def write_setup_sh(path, cfg, abs_rep_path):
    p, s = cfg["project_settings"], cfg["slurm_settings"]
    module_block = "\n".join(p.get("gmx_modules", []))
    analysis_script = os.path.join(os.getcwd(), p['scripts_dir'], p['density_analysis_script'])
    
    script_content = f"""#!/bin/bash
#SBATCH --job-name=S_{p['system_name']}
#SBATCH --account={s['account']}
#SBATCH --partition={s['partition_setup']}
#SBATCH --qos={s['qos_setup']}
#SBATCH --time={s['time_setup']}
#SBATCH --ntasks={s['ntasks_setup']}
#SBATCH --output={abs_rep_path}/setup.out
#SBATCH --error={abs_rep_path}/setup.err

{module_block}
source ~/miniconda3/bin/activate {p['conda_env']}

BASE="{abs_rep_path}"
GMX="gmx_mpi"

# 1. MIN
cd $BASE/1_min
if [ ! -f "min_out.gro" ]; then
    mpirun -np 1 $GMX grompp -f min.mdp -c start.gro -p topol.top -o min.tpr > $BASE/1_min/grompp_min.log 2>&1
    srun $GMX mdrun -deffnm min -c min_out.gro
else
    echo "Minimization already completed. Skipping."
fi

# 2. PRESS
cd $BASE/2_press
if [ ! -f "press_out.gro" ]; then
    mpirun -np 1 $GMX grompp -f press.mdp -c ../1_min/min_out.gro -p topol.top -o press.tpr > $BASE/2_press/grompp_press.log 2>&1
    srun $GMX mdrun -deffnm press -c press_out.gro -rdd 1.2
else
    echo "Pressure equilibration already completed. Skipping."
fi

# 3. ANNEAL
cd $BASE/3_anneal
if [ ! -f "anneal_out.gro" ]; then
    mpirun -np 1 $GMX grompp -f anneal.mdp -c ../2_press/press_out.gro -p topol.top -o anneal.tpr > $BASE/3_anneal/grompp_anneal.log 2>&1
    srun $GMX mdrun -deffnm anneal -c anneal_out.gro
else
    echo "Annealing already completed. Skipping."
fi

# 4. ANALYSIS
cd $BASE/3_anneal
if [ ! -f "$BASE/4_prod/start.gro" ]; then
    python3 {analysis_script} anneal.edr --num_replicas 1 --time_crop {p['density_analysis_time_crop']}
    [ -f "start_replica_1.gro" ] && mv start_replica_1.gro $BASE/4_prod/start.gro
else
    echo "Analysis and start.gro for production already exist. Skipping."
fi
"""
    with open(os.path.join(path, "run_setup.sh"), "w") as f: f.write(script_content)

def write_prod_sh(path, cfg, abs_rep_path, chunk_idx):
    p, s = cfg["project_settings"], cfg["slurm_settings"]
    module_block = "\n".join(p.get("gmx_modules", []))
    prod_path = os.path.join(abs_rep_path, "4_prod")
    script_content = f"""#!/bin/bash
#SBATCH --job-name=P{chunk_idx}_{p['system_name']}
#SBATCH --account={s['account']}
#SBATCH --partition={s['partition_prod']}
#SBATCH --qos={s['qos_prod']}
#SBATCH --time={s['time_prod']}
#SBATCH --ntasks={s['ntasks_prod']}
#SBATCH --output={prod_path}/chunk_{chunk_idx}.out
#SBATCH --error={prod_path}/chunk_{chunk_idx}.err

{module_block}
cd {prod_path}
if [ -f "nvt_run.cpt" ]; then
    srun gmx_mpi mdrun -v -deffnm nvt_run -cpi nvt_run.cpt -append -maxh 71.7
else
    mpirun -np 1 gmx_mpi grompp -f nvt_run.mdp -c start.gro -p topol.top -o nvt_run.tpr > {prod_path}/grompp_prod.log 2>&1
    srun gmx_mpi mdrun -v -deffnm nvt_run -maxh 71.7
fi
"""
    with open(os.path.join(path, f"run_prod_{chunk_idx}.sh"), "w") as f: f.write(script_content)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", nargs='?', default="config.toml")
    parser.add_argument("--dry-run", action="store_true", help="Print computed species counts and exit without running GROMACS or submitting jobs")
    args = parser.parse_args()

    if not os.path.exists(args.config): create_default_config(args.config)
    with open(args.config, "rb") as f: cfg = tomllib.load(f)

    scr = cfg["screening"]
    sp = normalize_species_config(cfg["species"])
    def to_list(val): return val if isinstance(val, list) else [val]
    group_defs, group_keys, active_species_keys = discover_groups(cfg, sp)
    order = active_species_keys

    # --- 1. DEFINE SCREENING AXES ---
    component_ratios = normalize_component_ratios(scr, group_keys)
    force_field_sets = normalize_force_field_sets(scr, order)

    # --- 2. CROSS COMBINATIONS ---
    # Cross (Component Ratio) x (Temperature) x (Force-Field Set)
    master_combos = list(itertools.product(component_ratios, to_list(scr["target_temps"]), force_field_sets))

    print(f"Setting {len(master_combos)} total variations...")

    if not args.dry_run:
        # --- GMX DETECTION ---
        gmx_path = cfg["project_settings"].get("gmx_executable_path")
        if not gmx_path or not os.path.exists(gmx_path): gmx_path = shutil.which("gmx_mpi")
        if not gmx_path: print("ERROR: GMX not found"); sys.exit(1)
        

    for ratio_entry, temp, ff_entry in master_combos:
        ratio_name = ratio_entry["name"]
        component_ratio = ratio_entry["weights"]
        ff_name = ff_entry["name"]
        active_itps = ff_entry["species_itps"]

        label = build_case_label(temp, ratio_name, ff_name, component_ratio, active_itps, order, group_keys)
        sys_path = os.path.abspath(f"data/{label}")
        
        counts, sizing_info = resolve_system_size(
            cfg, sp, order, group_defs, group_keys,
            component_ratio
        )

        if args.dry_run:
            print_dry_run_case(label, temp, ratio_name, ff_name, counts, order, group_keys, component_ratio, active_itps, sizing_info)
            continue

        for r in range(1, cfg['project_settings']['num_replicas'] + 1):
            rep_root = os.path.join(sys_path, f"rep_{r}")
            for stage in ["1_min", "2_press", "3_anneal", "4_prod"]:
                dest = os.path.join(rep_root, stage); os.makedirs(dest, exist_ok=True)
                for f in cfg["project_settings"]["global_ff_files"]: shutil.copy(f"inputs/{f}", f"{dest}/{f}")
                for k in order:
                    shutil.copy(f"inputs/{active_itps[k]}", f"{dest}/{active_itps[k]}")
                    shutil.copy(f"inputs/{sp[k]['gro']}", f"{dest}/{sp[k]['gro']}")
	    	    
                with open(os.path.join(dest, "topol.top"), "w") as top:
                    top.write('#include "forcefield.itp"\n')
                    for k in order: top.write(f'#include "{active_itps[k]}"\n')
                    top.write(f"\n[ system ]\n{label}\n\n[ molecules ]\n")
                    for k in order: top.write(f"{sp[k]['resname']:<15} {counts[k]}\n")
	    	
                tdir, sim = cfg["project_settings"]["template_dir"], cfg["simulation_settings"]
                if stage == "1_min": shutil.copy(f"{tdir}/min.mdp", f"{dest}/min.mdp")
                elif stage == "2_press":
                    modify_mdp(f"{tdir}/press.mdp", f"{dest}/press.mdp", {"ref_t": sim["press"]["ref_t"], "nsteps": sim["press"]["nsteps"], "dt": sim["press"]["dt"], "ref_p": sim["press"]["ref_p"]})
                elif stage == "3_anneal":
                    modify_mdp(f"{tdir}/anneal.mdp", f"{dest}/anneal.mdp", {
                        "ref_t": temp, "ref_p": sim["anneal"]["ref_p"], "nsteps": sim["anneal"]["nsteps"], "dt": sim["anneal"]["dt"],
                        "annealing-time": sim["anneal"]["annealing_time"], "annealing-temp": f"{sim['anneal']['temp_high']} {sim['anneal']['temp_high']} {temp} {temp}"
                    })
                elif stage == "4_prod":
                    modify_mdp(f"{tdir}/nvt_run.mdp", f"{dest}/nvt_run.mdp", {"ref_t": temp, "nsteps": sim["prod"]["nsteps"], "dt": sim["prod"]["dt"]})
	    	
            # --- BUILD BOX ---
            min_dir = os.path.join(rep_root, "1_min")
            ins_log_abs = os.path.join(min_dir, "insert-molecules.log")
            curr, box = None, sizing_info["box_size_nm"]
            with open(ins_log_abs, "w") as l: l.write("=== BOX LOG ===\n")

            for k in order:
                if counts[k] <= 0:
                    with open(ins_log_abs, "a") as l: l.write(f"\n--- Skipping {k} (nmol=0) ---\n")
                    continue
                out = f"tmp_{k}_{r}.gro"
                cmd = f"mpirun -np 1 {gmx_path} insert-molecules {'-box '+str(box)+' '+str(box)+' '+str(box) if curr is None else '-f '+curr} -ci inputs/{sp[k]['gro']} -nmol {counts[k]} -o {out} -seed {r*123}"
                res = subprocess.run(cmd.split(), capture_output=True, text=True)
                with open(ins_log_abs, "a") as l: l.write(f"\n--- Adding {k} ---\n" + res.stdout + res.stderr)
                if not os.path.exists(out): print(f"FAILED build for {label}. Check {ins_log_abs}"); sys.exit(1)
                curr = out

            if curr is None:
                print(f"FAILED build for {label}: no molecules were inserted. Check {ins_log_abs}")
                sys.exit(1)
            shutil.move(curr, os.path.join(min_dir, "start.gro"))
            for f in os.listdir("."): 
                if f.startswith("tmp_") and f.endswith(".gro"): os.remove(f)

            # --- SUBMIT SLURM ---
            write_setup_sh(rep_root, cfg, rep_root)
            sid = subprocess.run(["sbatch", "--parsable", os.path.join(rep_root, "run_setup.sh")], capture_output=True, text=True).stdout.strip()
            print(f"Launched {label} R{r}: Setup {sid}")
            prev = sid
            for c in range(1, cfg["project_settings"]["num_prod_chunks"] + 1):
                write_prod_sh(rep_root, cfg, rep_root, c)
                prev = subprocess.run(["sbatch", "--parsable", f"--dependency=afterany:{prev}", os.path.join(rep_root, f"run_prod_{c}.sh")], capture_output=True, text=True).stdout.strip()

if __name__ == "__main__": main()
