"""Molarity analysis using group counts and replica-specific box volumes."""

import os
import statistics

from ..geometry import parse_gro_box_volume_nm3

AVOGADRO = 6.02214076e23

DETAIL_HEADERS = [
    "System",
    "Replica",
    "Group",
    "Group Count",
    "Volume (nm^3)",
    "Volume (L)",
    "Molarity (mol/L)",
    "Status",
    "Structure File",
    "Error",
]

SUMMARY_HEADERS = [
    "System",
    "Group",
    "Configured Replicas",
    "Valid Replicas",
    "Mean Molarity (mol/L)",
    "Standard Deviation (mol/L)",
    "Standard Error (mol/L)",
    "Minimum (mol/L)",
    "Maximum (mol/L)",
    "Status",
]


def compute_molarity_mol_l(group_count, volume_nm3):
    if group_count < 0:
        raise ValueError("Group count must not be negative")
    if volume_nm3 <= 0:
        raise ValueError("Volume must be positive")
    return group_count / (AVOGADRO * volume_nm3 * 1e-24)


def _selected_groups(raw_groups, group_counts):
    if raw_groups in (None, "all"):
        return list(group_counts)
    if not isinstance(raw_groups, list) or not raw_groups:
        raise ValueError("analysis.molarity.groups must be 'all' or a non-empty list")
    unknown = [group for group in raw_groups if group not in group_counts]
    if unknown:
        raise KeyError(f"Unknown molarity analysis groups: {', '.join(unknown)}")
    return raw_groups


def run_molarity_analysis(case, settings):
    structure = str(settings.get("structure", "4_prod/start.gro"))
    groups = _selected_groups(settings.get("groups", "all"), case.group_counts)
    details = []

    for replica in range(1, case.configured_replicas + 1):
        structure_path = os.path.join(case.case_path, f"rep_{replica}", structure)
        relative_path = os.path.relpath(structure_path, case.case_path)
        try:
            volume_nm3 = parse_gro_box_volume_nm3(structure_path)
            error = ""
        except (OSError, ValueError) as exc:
            volume_nm3 = None
            error = str(exc)

        for group in groups:
            group_count = int(case.group_counts[group])
            molarity = (
                compute_molarity_mol_l(group_count, volume_nm3)
                if volume_nm3 is not None
                else None
            )
            details.append({
                "System": case.system,
                "Replica": replica,
                "Group": group,
                "Group Count": group_count,
                "Volume (nm^3)": volume_nm3,
                "Volume (L)": volume_nm3 * 1e-24 if volume_nm3 is not None else None,
                "Molarity (mol/L)": molarity,
                "Status": "OK" if molarity is not None else "Error",
                "Structure File": relative_path,
                "Error": error,
            })

    summaries = []
    for group in groups:
        values = [
            row["Molarity (mol/L)"]
            for row in details
            if row["Group"] == group and row["Molarity (mol/L)"] is not None
        ]
        valid_count = len(values)
        std = statistics.stdev(values) if valid_count > 1 else 0.0 if valid_count == 1 else None
        sem = std / valid_count**0.5 if std is not None and valid_count else None
        summaries.append({
            "System": case.system,
            "Group": group,
            "Configured Replicas": case.configured_replicas,
            "Valid Replicas": valid_count,
            "Mean Molarity (mol/L)": statistics.fmean(values) if values else None,
            "Standard Deviation (mol/L)": std,
            "Standard Error (mol/L)": sem,
            "Minimum (mol/L)": min(values) if values else None,
            "Maximum (mol/L)": max(values) if values else None,
            "Status": "Complete" if valid_count == case.configured_replicas else "Incomplete",
        })
    return DETAIL_HEADERS, details, SUMMARY_HEADERS, summaries
