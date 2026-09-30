#!/usr/bin/env python3
"""
Cell Scatter Analysis Pipeline
Main CLI entry point

Usage:
  python run_pipeline.py init --data <path>
  python run_pipeline.py config
  python run_pipeline.py run --all
  python run_pipeline.py status
  python run_pipeline.py results
"""

import click
from pipeline.cli import (
    init_experiment,
    configure_experiment,
    run_stages,
    show_status,
    show_results,
    clean_outputs,
    validate_experiment
)


@click.group()
@click.version_option(version='1.0.0')
def cli():
    """Cell Scatter Analysis Pipeline

    A modular 7-stage pipeline for analyzing microscopy images to quantify
    cellular response to treatments by measuring particle scattering patterns.
    """
    pass


@cli.command()
@click.option('--data', required=True, type=click.Path(exists=True),
              help='Path to data directory')
@click.option('--name', default='Experiment', help='Experiment name')
@click.option('--description', default='', help='Experiment description')
@click.option('--output', default='outputs', help='Output directory')
@click.option('--config', 'config_file', type=click.Path(exists=True),
              help='Use existing config file')
@click.option('--template', type=click.Choice(['temporal', 'marker', 'custom']),
              help='Config template')
def init(data, name, description, output, config_file, template):
    """Initialize new experiment"""
    init_experiment(data, name, description, output, config_file, template)


@cli.command()
@click.option('--interactive', is_flag=True, default=True,
              help='Interactive configuration (default)')
@click.option('--edit', is_flag=True, help='Edit config file')
@click.option('--show', is_flag=True, help='Show current config')
@click.option('--validate', 'validate_only', is_flag=True, help='Validate config')
def config(interactive, edit, show, validate_only):
    """Configure experiment channels and comparisons"""
    configure_experiment(interactive, edit, show, validate_only)


@cli.command()
@click.option('--stage', type=int, help='Run specific stage (1-7)')
@click.option('--all', 'run_all', is_flag=True, help='Run all stages')
@click.option('--from', 'from_stage', type=int, help='Run from stage N')
@click.option('--to', 'to_stage', type=int, help='Run to stage N')
@click.option('--force', is_flag=True, help='Force rerun (ignore cache)')
@click.option('--dry-run', is_flag=True, help='Dry run (show plan)')
def run(stage, run_all, from_stage, to_stage, force, dry_run):
    """Run pipeline stages"""
    run_stages(stage, run_all, from_stage, to_stage, force, dry_run)


@cli.command()
@click.option('--detailed', is_flag=True, help='Detailed status')
@click.option('--json', 'as_json', is_flag=True, help='JSON output')
def status(detailed, as_json):
    """Show pipeline status"""
    show_status(detailed, as_json)


@cli.command()
@click.option('--summary', is_flag=True, help='Show summary')
@click.option('--cell', help='Show specific cell results')
@click.option('--comparison', help='Show specific comparison results')
@click.option('--export', type=click.Choice(['csv', 'json', 'excel']),
              help='Export results')
@click.option('--report', is_flag=True, help='Generate HTML report')
def results(summary, cell, comparison, export, report):
    """View and export results"""
    show_results(summary, cell, comparison, export, report)


@cli.command()
@click.option('--stage', type=int, help='Clean specific stage')
@click.option('--all', 'clean_all', is_flag=True, help='Clean all outputs')
@click.option('--cache', is_flag=True, help='Clean cache only')
@click.option('--confirm', is_flag=True, help='Skip confirmation')
def clean(stage, clean_all, cache, confirm):
    """Clean outputs and reset pipeline"""
    clean_outputs(stage, clean_all, cache, confirm)


@cli.command()
@click.option('--data', 'data_only', is_flag=True, help='Validate data')
@click.option('--config', 'config_only', is_flag=True, help='Validate config')
@click.option('--outputs', is_flag=True, help='Validate outputs')
@click.option('--all', 'validate_all', is_flag=True, help='Validate everything')
def validate(data_only, config_only, outputs, validate_all):
    """Validate experiment setup and data"""
    validate_experiment(data_only, config_only, outputs, validate_all)


if __name__ == '__main__':
    cli()
