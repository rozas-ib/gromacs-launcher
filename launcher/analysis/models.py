"""Common result models for post-simulation analyses."""

from dataclasses import dataclass


@dataclass(frozen=True)
class AnalysisCase:
    system: str
    case_path: str
    configured_replicas: int
    group_counts: dict
