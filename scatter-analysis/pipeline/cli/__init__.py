"""
Command Line Interface

CLI commands for the pipeline:
- init: Initialize experiment
- config: Configure channels and comparisons
- run: Execute pipeline stages
- status: Check pipeline status
- results: View and export results
- validate: Validate experiment setup
- clean: Clean outputs
"""

from .commands import (
    init_experiment,
    configure_experiment,
    run_stages,
    show_status,
    show_results,
    clean_outputs,
    validate_experiment
)

__all__ = [
    'init_experiment',
    'configure_experiment',
    'run_stages',
    'show_status',
    'show_results',
    'clean_outputs',
    'validate_experiment'
]
