import itertools
import math
import os
from dataclasses import dataclass

from .paths import get_case_dir


@dataclass(frozen=True)
class CaseContext:
    label: str
    sys_path: str
    job_log_path: str
    ratio_entry: dict
    temp: float
    ff_entry: dict
    active_itps: dict


def to_list(value):
    return value if isinstance(value, list) else [value]


def read_gro_species_metadata(gro_path):
    with open(gro_path, "r", encoding="utf-8") as f:
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
        for key in species_keys:
            value = entry[key]
            if not isinstance(value, str) or not value.strip():
                raise ValueError(
                    f"screening.force_field '{name}' entry for species '{key}' must be a non-empty string"
                )
            species_itps[key] = value.strip()
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
        for key in rank[:remainder]:
            counts[key] += 1
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
            raise ValueError(
                "system_sizing.atom_number_density_atoms_per_nm3 must be > 0 when estimate_box_from_atoms = true"
            )
        volume_nm3 = float(total_atoms) / atom_density
        scale_factor = float(project_cfg.get("box_scale_factor", 1.0))
        return (scale_factor * volume_nm3) ** (1.0 / 3.0), "estimated_and_scaled_from_atom_density"
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
    itp_label = "_".join(active_itps[k].split(".")[0] for k in species_order)
    return f"T{fmt_num(temp)}_{ratio_name}_{ff_name}_{comp_mix_label}_{itp_label}"


def print_dry_run_case(label, temp, ratio_name, ff_name, counts, order, group_keys, component_ratio, active_itps, sizing_info):
    total_system = sum(counts[k] for k in order)
    component_weights = ", ".join(f"{k}={fmt_num(component_ratio[k])}" for k in group_keys)
    species_counts = ", ".join(f"{k}={counts[k]}" for k in order)
    print(f"{'=' * 100}\n[DRY-RUN] {label}\n{'=' * 100}")
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


def format_final_molar_ratio(component_counts, group_keys):
    values = [int(component_counts.get(k, 0)) for k in group_keys]
    positive = [v for v in values if v > 0]
    if not positive:
        return "all components are zero"
    gcd_val = positive[0]
    for value in positive[1:]:
        gcd_val = math.gcd(gcd_val, value)
    reduced = {k: (int(component_counts.get(k, 0)) // gcd_val) for k in group_keys}
    total = sum(values)
    fractions = {k: (int(component_counts.get(k, 0)) / total if total > 0 else 0.0) for k in group_keys}
    reduced_str = ":".join(str(reduced[k]) for k in group_keys)
    frac_str = ", ".join(f"{k}={fractions[k]:.6f}" for k in group_keys)
    return f"reduced={reduced_str} ({', '.join(group_keys)}), mole_fractions=({frac_str})"


def build_case_contexts(master_combos, order, group_keys, output_root):
    contexts = []
    for ratio_entry, temp, ff_entry in master_combos:
        ratio_name = ratio_entry["name"]
        component_ratio = ratio_entry["weights"]
        ff_name = ff_entry["name"]
        active_itps = ff_entry["species_itps"]
        label = build_case_label(temp, ratio_name, ff_name, component_ratio, active_itps, order, group_keys)
        sys_path = get_case_dir(output_root, label)
        contexts.append(
            CaseContext(
                label=label,
                sys_path=sys_path,
                job_log_path=os.path.join(sys_path, "launch.log"),
                ratio_entry=ratio_entry,
                temp=temp,
                ff_entry=ff_entry,
                active_itps=active_itps,
            )
        )
    return contexts


def build_master_combos(cfg, group_keys, species_order):
    scr = cfg["screening"]
    component_ratios = normalize_component_ratios(scr, group_keys)
    force_field_sets = normalize_force_field_sets(scr, species_order)
    return list(itertools.product(component_ratios, to_list(scr["target_temps"]), force_field_sets))

