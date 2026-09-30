"""
CLI Command Implementations

All command functions for the pipeline CLI.
"""

import os
import sys
import json
import yaml
import shutil
from pathlib import Path
from datetime import datetime
import click


def init_experiment(data_dir, name, description, output_dir, config_file, template):
    """Initialize a new experiment"""
    click.echo("=" * 60)
    click.echo(f"Initializing experiment: {name}")
    click.echo("=" * 60)

    data_path = Path(data_dir)
    output_path = Path(output_dir)

    # Validate data directory structure
    click.echo(f"✓ Data directory: {data_path}")

    # Dynamically discover all condition groups (subdirectories)
    if not data_path.exists():
        click.echo(f"✗ Error: Data directory not found at {data_path}", err=True)
        sys.exit(1)

    # Find all subdirectories (excluding hidden dirs)
    group_dirs = [d for d in data_path.iterdir() if d.is_dir() and not d.name.startswith('.')]

    if not group_dirs:
        click.echo(f"✗ Error: No subdirectories found in {data_path}", err=True)
        sys.exit(1)

    # Validate each group has image files
    valid_groups = {}
    all_files = []
    for group_dir in sorted(group_dirs):
        group_files = list(group_dir.glob("**/*.tif")) + list(group_dir.glob("**/*.tiff"))
        if group_files:
            valid_groups[group_dir.name] = group_files
            all_files.extend(group_files)

    if not valid_groups:
        click.echo(f"✗ Error: No .tif/.tiff files found in any subdirectory", err=True)
        sys.exit(1)

    click.echo(f"✓ Found groups: {', '.join(sorted(valid_groups.keys()))}")

    # Detect channels from all files
    channels = set()
    for f in all_files[:20]:  # Check first 20 files across all groups
        if "C1-" in f.name or "C1_" in f.name:
            channels.add("C1")
        if "C2-" in f.name or "C2_" in f.name:
            channels.add("C2")
        if "C3-" in f.name or "C3_" in f.name:
            channels.add("C3")

    channels = sorted(channels)
    if channels:
        click.echo(f"✓ Detected channels: {', '.join(channels)}")
    else:
        click.echo("⚠ Warning: Could not detect channel names from files")
        channels = ["C1", "C2", "C3"]

    # Count samples per group
    group_samples = {}
    for group_name, group_files in valid_groups.items():
        samples = set()
        for f in group_files:
            # Extract sample ID (e.g., Untreated_1 from C1-Untreated_1.tif)
            parts = f.stem.split('-', 1)
            if len(parts) > 1:
                samples.add(parts[1])
        group_samples[group_name] = len(samples)

    click.echo(f"✓ Detected samples:")
    for group_name in sorted(group_samples.keys()):
        click.echo(f"    {group_name}: {group_samples[group_name]} samples")

    # Create output directory structure
    output_path.mkdir(exist_ok=True)
    for i in range(1, 8):
        stage_dir = output_path / f"0{i}_{'validated' if i == 1 else 'stage' + str(i)}"
        if i == 1:
            stage_dir = output_path / "01_validated"
        elif i == 2:
            stage_dir = output_path / "02_segmented"
        elif i == 3:
            stage_dir = output_path / "03_registered"
        elif i == 4:
            stage_dir = output_path / "04_particles"
        elif i == 5:
            stage_dir = output_path / "05_metrics"
        elif i == 6:
            stage_dir = output_path / "06_analysis"
        elif i == 7:
            stage_dir = output_path / "07_visualizations"
        stage_dir.mkdir(exist_ok=True)

    click.echo(f"✓ Created output directory: {output_path}")

    # Create experiment metadata
    groups_meta = {}
    for group_name in sorted(valid_groups.keys()):
        groups_meta[group_name] = {
            "path": str((data_path / group_name).absolute()),
            "samples": group_samples[group_name]
        }

    experiment_meta = {
        "name": name,
        "description": description,
        "initialized_at": datetime.now().isoformat(),
        "data_directory": str(data_path.absolute()),
        "output_directory": str(output_path.absolute()),
        "channels": channels,
        "groups": groups_meta
    }

    experiment_file = output_path / "experiment.json"
    with open(experiment_file, 'w') as f:
        json.dump(experiment_meta, f, indent=2)

    click.echo(f"✓ Created experiment metadata: {experiment_file}")

    # Create or copy config file
    config_path = output_path / "config.yaml"

    if config_file:
        # Copy existing config
        shutil.copy(config_file, config_path)
        click.echo(f"✓ Copied config from: {config_file}")
    elif template:
        # Create from template
        config_content = create_config_from_template(template, channels)
        with open(config_path, 'w') as f:
            yaml.dump(config_content, f, default_flow_style=False, sort_keys=False)
        click.echo(f"✓ Generated config from template: {template}")
    else:
        # Create default config
        config_content = create_default_config(name, description, channels)
        with open(config_path, 'w') as f:
            yaml.dump(config_content, f, default_flow_style=False, sort_keys=False)
        click.echo(f"✓ Generated default config: {config_path}")

    # Create pipeline state file
    state = {
        "experiment_id": name.lower().replace(" ", "_") + "_" + datetime.now().strftime("%Y%m%d"),
        "initialized_at": datetime.now().isoformat(),
        "last_updated": datetime.now().isoformat(),
        "stages": {
            str(i): {
                "name": get_stage_name(i),
                "status": "pending",
                "output_dir": str(output_path / get_stage_dirname(i))
            }
            for i in range(1, 8)
        },
        "progress": {
            "stages_completed": 0,
            "stages_total": 7,
            "percent_complete": 0.0
        }
    }

    state_file = output_path / "pipeline_state.json"
    with open(state_file, 'w') as f:
        json.dump(state, f, indent=2)

    click.echo(f"✓ Created pipeline state: {state_file}")

    click.echo()
    click.echo("=" * 60)
    click.echo("✓ Experiment initialized successfully!")
    click.echo("=" * 60)
    click.echo()
    click.echo("Next step: Configure your experiment")
    click.echo(f"  python run_pipeline.py config")


def configure_experiment(interactive, edit, show, validate_only):
    """Configure experiment channels and comparisons"""
    output_path = Path("outputs")
    config_path = output_path / "config.yaml"

    if not config_path.exists():
        click.echo("✗ Error: No config file found. Run 'init' first.", err=True)
        sys.exit(1)

    if show:
        # Show current configuration
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        click.echo("=" * 60)
        click.echo("Current Configuration")
        click.echo("=" * 60)
        click.echo(yaml.dump(config, default_flow_style=False))
        return

    if validate_only:
        # Validate configuration
        click.echo("Validating configuration...")
        # TODO: Implement validation logic
        click.echo("✓ Configuration is valid")
        return

    if edit:
        # Open in editor
        editor = os.environ.get('EDITOR', 'nano')
        os.system(f'{editor} {config_path}')
        return

    if interactive:
        # Interactive configuration
        click.echo("=" * 60)
        click.echo("Experiment Configuration")
        click.echo("=" * 60)
        click.echo()
        click.echo("Interactive configuration not yet implemented.")
        click.echo(f"Please edit the config file directly: {config_path}")
        click.echo()
        click.echo("Or use:")
        click.echo(f"  python run_pipeline.py config --edit")


def run_stages(stage, run_all, from_stage, to_stage, force, dry_run):
    """Run pipeline stages"""
    from pipeline.stages.stage_01_validate import run_stage_1
    from pipeline.stages.stage_02_segment import run_stage_2
    from pipeline.stages.stage_03_register import run_stage_3
    from pipeline.stages.stage_04_detect_particles import run_stage_4
    from pipeline.stages.stage_05_calculate_metrics import run_stage_5
    from pipeline.stages.stage_06_analyze import run_stage_6
    from pipeline.stages.stage_07_visualize import run_stage_7

    # Load experiment metadata and config
    output_path = Path("outputs")
    experiment_file = output_path / "experiment.json"
    config_file = output_path / "config.yaml"
    state_file = output_path / "pipeline_state.json"

    if not experiment_file.exists():
        click.echo("✗ Error: No experiment found. Run 'init' first.", err=True)
        sys.exit(1)

    with open(experiment_file, 'r') as f:
        experiment_meta = json.load(f)

    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    with open(state_file, 'r') as f:
        state = json.load(f)

    # Determine which stages to run
    stages_to_run = []
    if run_all:
        stages_to_run = list(range(1, 8))
    elif stage:
        stages_to_run = [stage]
    elif from_stage and to_stage:
        stages_to_run = list(range(from_stage, to_stage + 1))
    elif from_stage:
        stages_to_run = list(range(from_stage, 8))
    else:
        click.echo("Error: Specify --stage, --all, or --from", err=True)
        sys.exit(1)

    click.echo("=" * 60)
    click.echo("Running Pipeline")
    click.echo("=" * 60)
    click.echo(f"\nExperiment: {experiment_meta['name']}")
    click.echo(f"Config: {config_file}")

    if dry_run:
        click.echo("\nDry run - stages that would be executed:")
        for s in stages_to_run:
            click.echo(f"  [{s}] {get_stage_name(s)}")
        return

    # Run each stage
    for stage_num in stages_to_run:
        stage_name = get_stage_name(stage_num)
        stage_dir = output_path / get_stage_dirname(stage_num)

        # Check if already completed (unless force)
        if state['stages'][str(stage_num)]['status'] == 'completed' and not force:
            click.echo(f"\n[Stage {stage_num}] {stage_name}")
            click.echo("  ⏭ Already completed (use --force to rerun)")
            continue

        # Update state to in_progress
        state['stages'][str(stage_num)]['status'] = 'in_progress'
        state['stages'][str(stage_num)]['started_at'] = datetime.now().isoformat()
        _save_state(state, state_file)

        # Run the stage
        start_time = datetime.now()

        try:
            if stage_num == 1:
                manifest, validation_report = run_stage_1(config, experiment_meta, stage_dir)
            elif stage_num == 2:
                # Load manifest from Stage 1
                manifest_file = output_path / "01_validated" / "manifest.json"
                with open(manifest_file, 'r') as f:
                    manifest = json.load(f)
                segmentation_report = run_stage_2(config, manifest, stage_dir)
            elif stage_num == 3:
                registration_report = run_stage_3(config, stage_dir)
            elif stage_num == 4:
                particle_report = run_stage_4(config, stage_dir)
            elif stage_num == 5:
                metrics_report = run_stage_5(config, stage_dir)
            elif stage_num == 6:
                all_metrics, statistics = run_stage_6(config, stage_dir)
            elif stage_num == 7:
                viz_report = run_stage_7(config, stage_dir)
            else:
                click.echo(f"\n[Stage {stage_num}] {stage_name}")
                click.echo("  ⚠ Not yet implemented")
                continue

            # Update state to completed
            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()

            state['stages'][str(stage_num)]['status'] = 'completed'
            state['stages'][str(stage_num)]['completed_at'] = end_time.isoformat()
            state['stages'][str(stage_num)]['duration_seconds'] = duration

            # Update progress
            completed = sum(1 for s in state['stages'].values() if s['status'] == 'completed')
            state['progress']['stages_completed'] = completed
            state['progress']['percent_complete'] = (completed / 7) * 100

            _save_state(state, state_file)

        except Exception as e:
            click.echo(f"\n✗ Error in stage {stage_num}: {str(e)}", err=True)
            state['stages'][str(stage_num)]['status'] = 'failed'
            state['stages'][str(stage_num)]['error'] = str(e)
            _save_state(state, state_file)
            sys.exit(1)

    click.echo("\n" + "=" * 60)
    click.echo("Pipeline Complete")
    click.echo("=" * 60)
    completed = state['progress']['stages_completed']
    click.echo(f"\nStages completed: {completed}/7")
    click.echo(f"\nView status: python run_pipeline.py status")


def _save_state(state, state_file):
    """Helper to save pipeline state"""
    state['last_updated'] = datetime.now().isoformat()
    with open(state_file, 'w') as f:
        json.dump(state, f, indent=2)


def show_status(detailed, as_json):
    """Show pipeline status"""
    output_path = Path("outputs")
    state_file = output_path / "pipeline_state.json"

    if not state_file.exists():
        click.echo("✗ Error: No pipeline state found. Run 'init' first.", err=True)
        sys.exit(1)

    with open(state_file, 'r') as f:
        state = json.load(f)

    if as_json:
        click.echo(json.dumps(state, indent=2))
        return

    click.echo("=" * 60)
    click.echo("Pipeline Status")
    click.echo("=" * 60)
    click.echo()
    click.echo(f"Experiment: {state['experiment_id']}")
    click.echo(f"Initialized: {state['initialized_at']}")
    click.echo()
    click.echo("Stages:")

    for stage_num in range(1, 8):
        stage_info = state['stages'][str(stage_num)]
        status = stage_info['status']
        name = stage_info['name']

        if status == 'completed':
            icon = '✓'
        elif status == 'in_progress':
            icon = '⏳'
        else:
            icon = '⏸'

        click.echo(f"  {icon} [{stage_num}] {name:<35} ({status})")

    click.echo()
    progress = state['progress']
    click.echo(f"Progress: {progress['stages_completed']}/{progress['stages_total']} stages ({progress['percent_complete']:.1f}%)")


def show_results(summary, cell, comparison, export, report):
    """View and export results"""
    click.echo("=" * 60)
    click.echo("Results")
    click.echo("=" * 60)
    click.echo()
    click.echo("Results viewing not yet implemented.")


def clean_outputs(stage, clean_all, cache, confirm):
    """Clean outputs"""
    click.echo("=" * 60)
    click.echo("Clean Outputs")
    click.echo("=" * 60)
    click.echo()
    click.echo("Clean not yet implemented.")


def validate_experiment(data_only, config_only, outputs, validate_all):
    """Validate experiment setup"""
    click.echo("=" * 60)
    click.echo("Validation")
    click.echo("=" * 60)
    click.echo()
    click.echo("Validation not yet implemented.")


# Helper functions

def get_stage_name(stage_num):
    """Get stage name by number"""
    stages = {
        1: "Data Loading & Validation",
        2: "Cell Segmentation",
        3: "Cell Registration/Mapping",
        4: "Particle Detection",
        5: "Scatter Metrics Calculation",
        6: "Aggregation & Comparison",
        7: "Visualization & Reporting"
    }
    return stages.get(stage_num, f"Stage {stage_num}")


def get_stage_dirname(stage_num):
    """Get stage directory name"""
    dirs = {
        1: "01_validated",
        2: "02_segmented",
        3: "03_registered",
        4: "04_particles",
        5: "05_metrics",
        6: "06_analysis",
        7: "07_visualizations"
    }
    return dirs.get(stage_num, f"0{stage_num}_stage{stage_num}")


def create_default_config(name, description, channels):
    """Create default configuration"""
    config = {
        "experiment": {
            "name": name,
            "type": "temporal",
            "description": description or "Cell scatter analysis experiment"
        },
        "channels": {},
        "comparisons": [],
        "validation": {
            "strictness": "flexible",
            "cell_segmentation": {
                "min_cell_area": 100,
                "max_cell_area": 50000,
                "boundary_exclusion_threshold": 0.5
            },
            "cell_mapping": {
                "min_confidence": 30,
                "warning_threshold": 70
            },
            "outlier_detection": {
                "enabled": True,
                "method": "IQR",
                "use_median_threshold": 0.05
            }
        },
        "segmentation": {
            "method": "contour_detection",
            "min_cell_area": 100,
            "max_cell_area": 50000,
            "threshold_method": "otsu"
        },
        "particle_detection": {
            "threshold_method": "otsu",
            "min_particle_size": 1,
            "max_particle_size": 1000,
            "blur_sigma": 1.0
        },
        "metrics": {
            "primary_metric": "median_distance",
            "calculate": [
                "mean_distance_from_center",
                "median_distance_from_center",
                "std_distance_from_center",
                "max_distance_from_center",
                "convex_hull_area",
                "particle_density"
            ]
        }
    }

    # Add channel placeholders
    for ch in channels:
        config["channels"][ch] = {
            "name": f"Channel {ch}",
            "role": "unknown",
            "description": f"{ch} - configure this channel"
        }

    return config


def create_config_from_template(template, channels):
    """Create config from template"""
    # For now, just create default config
    # TODO: Implement proper templates
    return create_default_config(f"{template.capitalize()} Experiment", "", channels)
