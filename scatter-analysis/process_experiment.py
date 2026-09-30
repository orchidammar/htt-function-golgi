#!/usr/bin/env python3
"""
Process a single experiment through the pipeline (stages 1-5)
Usage: python process_experiment.py <experiment_dir> <output_dir>
"""

import sys
import yaml
from pathlib import Path
from pipeline.stages.stage_01_validate import run_stage_1
from pipeline.stages.stage_02_segment import run_stage_2
from pipeline.stages.stage_03_register import run_stage_3
from pipeline.stages.stage_04_detect_particles import run_stage_4
from pipeline.stages.stage_05_calculate_metrics import run_stage_5

def create_config(experiment_dir, output_dir):
    """Create config for an experiment"""
    config = {
        'experiment': {
            'name': Path(experiment_dir).name,
            'data_directory': str(experiment_dir)
        },
        'channels': {
            'c1': {
                'name': 'Structure Y',
                'color': 'red',
                'biological_role': 'Reference structure'
            },
            'c2': {
                'name': 'Protein X',
                'color': 'green',
                'biological_role': 'Target protein'
            },
            'c3': {
                'name': 'Positioning',
                'color': 'blue',
                'biological_role': 'Cell positioning marker'
            }
        },
        'groups': {
            'comparison_baseline': 'Controls'
        },
        'segmentation': {
            'min_cell_size': 100,
            'max_cell_size': 10000,
            'boundary_exclusion': 0.1
        },
        'particle_detection': {
            'min_particle_size': 3,
            'max_particle_size': 50,
            'intensity_threshold_method': 'otsu'
        }
    }

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    config_file = output_path / 'config.yaml'
    with open(config_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return config, output_path

def main():
    if len(sys.argv) != 3:
        print("Usage: python process_experiment.py <experiment_dir> <output_dir>")
        sys.exit(1)

    experiment_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])

    if not experiment_dir.exists():
        print(f"Error: {experiment_dir} does not exist")
        sys.exit(1)

    print(f"\n{'='*80}")
    print(f"Processing {experiment_dir.name}")
    print(f"Output: {output_dir}")
    print(f"{'='*80}\n")

    # Create config
    config, output_path = create_config(experiment_dir, output_dir)

    # Stage 1: Validate
    print("\n[Stage 1] Validating data...")
    stage_01_dir = output_path / "01_validated"
    experiment_meta = {
        'name': config['experiment']['name'],
        'data_directory': config['experiment']['data_directory'],
        'channels': ['C1', 'C2', 'C3']  # Channel prefixes in filenames
    }
    manifest, validation_report = run_stage_1(config, experiment_meta, stage_01_dir)

    # Stage 2: Segment
    print("\n[Stage 2] Segmenting cells...")
    stage_02_dir = output_path / "02_segmented"
    run_stage_2(config, manifest, stage_02_dir)

    # Stage 3: Register
    print("\n[Stage 3] Registering channels...")
    stage_03_dir = output_path / "03_registered"
    run_stage_3(config, stage_03_dir)

    # Stage 4: Detect particles
    print("\n[Stage 4] Detecting particles...")
    stage_04_dir = output_path / "04_particles"
    run_stage_4(config, stage_04_dir)

    # Stage 5: Calculate metrics (including scatter indices)
    print("\n[Stage 5] Calculating metrics...")
    stage_05_dir = output_path / "05_metrics"
    run_stage_5(config, stage_05_dir)

    print(f"\n{'='*80}")
    print(f"✅ {experiment_dir.name} processing complete!")
    print(f"Output: {output_dir}")
    print(f"{'='*80}\n")

if __name__ == '__main__':
    main()
