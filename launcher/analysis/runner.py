"""Configuration-driven dispatcher for post-simulation analyses."""

import os

from ..case_matrix import (
    build_case_contexts,
    build_master_combos,
    discover_groups,
    normalize_species_config,
    resolve_system_size,
)
from ..config import load_config
from ..paths import get_inputs_dir, get_output_root
from .models import AnalysisCase
from .registry import ANALYSIS_RUNNERS
from .reporting import write_csv


def _build_cases(cfg):
    species_cfg = normalize_species_config(cfg["species"], inputs_dir=get_inputs_dir(cfg))
    group_defs, group_keys, order = discover_groups(cfg, species_cfg)
    output_root = get_output_root(cfg)
    contexts = build_case_contexts(build_master_combos(cfg, group_keys, order), order, group_keys, output_root)
    replica_count = int(cfg["project_settings"]["num_replicas"])
    for context in contexts:
        _, sizing_info = resolve_system_size(
            cfg,
            species_cfg,
            order,
            group_defs,
            group_keys,
            context.ratio_entry["weights"],
            context.sizing_entry,
        )
        yield AnalysisCase(
            system=context.label,
            case_path=context.sys_path,
            configured_replicas=replica_count,
            group_counts=sizing_info["component_counts"],
        )


def run_analysis(config_path):
    cfg = load_config(config_path)
    analysis_cfg = cfg.get("analysis", {})
    if not analysis_cfg.get("enabled", False):
        raise ValueError("Post-simulation analysis is disabled; set [analysis].enabled = true")

    output_root = get_output_root(cfg)
    output_dir = os.path.join(output_root, str(analysis_cfg.get("output_subdir", "analysis")))
    enabled = [
        (name, settings)
        for name, settings in analysis_cfg.items()
        if name in ANALYSIS_RUNNERS and isinstance(settings, dict) and settings.get("enabled", False)
    ]
    if not enabled:
        raise ValueError("No post-simulation analyses are enabled")

    written = []
    cases = list(_build_cases(cfg))
    for name, settings in enabled:
        detail_headers = summary_headers = None
        details = []
        summaries = []
        for case in cases:
            detail_headers, case_details, summary_headers, case_summaries = ANALYSIS_RUNNERS[name](case, settings)
            details.extend(case_details)
            summaries.extend(case_summaries)
        detail_path = write_csv(os.path.join(output_dir, f"{name}_per_replica.csv"), detail_headers, details)
        summary_path = write_csv(os.path.join(output_dir, f"{name}_summary.csv"), summary_headers, summaries)
        written.extend([detail_path, summary_path])
    return written
