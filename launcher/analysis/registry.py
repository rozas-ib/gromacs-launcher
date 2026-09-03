"""Registry of available post-simulation analysis wrappers."""

from .analyses.molarity import run_molarity_analysis


ANALYSIS_RUNNERS = {
    "molarity": run_molarity_analysis,
}
