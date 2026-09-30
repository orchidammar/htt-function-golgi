"""
Stage 3: Cell Registration/Mapping

Purpose:
    Verify and document cell correspondence across channels.

SIMPLIFIED APPROACH:
    Since we use C1-as-reference strategy with SAME bounding boxes:
    - Cell 1 in C1 = Cell 1 in C2 = Cell 1 in C3 (same spatial region)
    - No complex matching needed
    - Cells are already aligned by design

Input:
    - outputs/02_segmented/[group]/[sample_id]/sample_metadata.json

Output:
    - outputs/03_registered/[group]/[sample_id]/cell_mapping.json
    - outputs/03_registered/[group]/[sample_id]/cell_[n]_metadata.json
"""

import json
import cv2
import numpy as np
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
import click


class CellRegistrar:
    """Registers cells across channels (simplified since already aligned)"""

    def __init__(self, config, output_dir):
        self.config = config
        self.output_dir = Path(output_dir)
        self.segmented_dir = Path("outputs/02_segmented")

        # Statistics
        self.stats = {
            'samples_processed': 0,
            'samples_failed': 0,
            'total_cells': 0,
            'cells_mapped': 0,
            'cells_with_issues': 0
        }

        self.registration_report = {
            'created_at': datetime.now().isoformat(),
            'samples': {},
            'errors': [],
            'warnings': []
        }

    def run(self):
        """Execute Stage 3: Cell Registration/Mapping"""
        click.echo("\n" + "=" * 60)
        click.echo("[Stage 3] Cell Registration/Mapping")
        click.echo("=" * 60)

        # Load segmentation report to get list of samples
        seg_report_path = self.segmented_dir / "segmentation_report.json"
        with open(seg_report_path, 'r') as f:
            seg_report = json.load(f)

        samples = seg_report.get('samples', {})

        click.echo(f"\n→ Registering cells across channels for {len(samples)} samples...")
        click.echo("   (Using C1-as-reference strategy - cells already aligned)")

        for sample_id, sample_data in tqdm(samples.items(), desc="Registering"):
            try:
                self._register_sample(sample_id, sample_data)
                self.stats['samples_processed'] += 1
            except Exception as e:
                click.echo(f"\n✗ Error registering {sample_id}: {str(e)}")
                self.stats['samples_failed'] += 1
                self.registration_report['errors'].append({
                    'sample': sample_id,
                    'error': str(e)
                })

        # Save registration report
        self._save_report()

        # Print results
        self._print_results()

        return self.registration_report

    def _register_sample(self, sample_id, sample_data):
        """Register cells for a single sample"""

        group = sample_data['group']
        sample_name = sample_data['sample_name']

        # Create output directory
        sample_output_dir = self.output_dir / group / sample_name
        sample_output_dir.mkdir(parents=True, exist_ok=True)

        cells = sample_data.get('cells', [])
        self.stats['total_cells'] += len(cells)

        # Create cell mapping
        cell_mapping = {
            'sample_id': sample_id,
            'group': group,
            'sample_name': sample_name,
            'total_cells': len(cells),
            'cells': {}
        }

        # For each cell, verify correspondence across channels
        for cell_data in cells:
            cell_id = cell_data['cell_id']
            channels_data = cell_data['channels']

            # Since we used same bbox for all channels, correspondence is trivial
            cell_info = {
                'cell_id': cell_id,
                'bbox': cell_data['bbox'],
                'center': cell_data['center'],
                'area': cell_data['area'],
                'correspondence': 'DIRECT',  # Same bbox across all channels
                'confidence': 100.0,  # 100% confidence - by design
                'status': 'ALIGNED',
                'channels': {}
            }

            # Verify each channel and perform sanity checks
            for channel, ch_data in channels_data.items():
                # Load cell image to verify
                cell_img_path = Path(ch_data['path'])

                if not cell_img_path.exists():
                    self.registration_report['errors'].append({
                        'code': 'E401',
                        'sample': sample_id,
                        'cell': cell_id,
                        'channel': channel,
                        'message': f'Cell image not found: {cell_img_path}'
                    })
                    cell_info['status'] = 'ERROR'
                    continue

                # Check if region has signal (not completely empty)
                mean_intensity = ch_data.get('mean_intensity', 0)

                if mean_intensity < 0.1:  # Very low signal
                    self.registration_report['warnings'].append({
                        'code': 'E402',
                        'sample': sample_id,
                        'cell': cell_id,
                        'channel': channel,
                        'message': f'Cell region appears empty (mean intensity: {mean_intensity:.2f})'
                    })
                    self.stats['cells_with_issues'] += 1

                # Add channel info
                cell_info['channels'][channel] = {
                    'path': str(cell_img_path),
                    'shape': ch_data['shape'],
                    'mean_intensity': mean_intensity,
                    'max_intensity': ch_data.get('max_intensity', 0),
                    'signal_quality': 'GOOD' if mean_intensity > 1.0 else 'LOW'
                }

            # Save individual cell metadata
            cell_metadata_path = sample_output_dir / f"{cell_id}_metadata.json"
            with open(cell_metadata_path, 'w') as f:
                json.dump(cell_info, f, indent=2)

            # Add to mapping
            cell_mapping['cells'][cell_id] = {
                'correspondence': cell_info['correspondence'],
                'confidence': cell_info['confidence'],
                'status': cell_info['status'],
                'metadata_path': str(cell_metadata_path)
            }

            if cell_info['status'] == 'ALIGNED':
                self.stats['cells_mapped'] += 1

        # Save cell mapping for this sample
        mapping_path = sample_output_dir / 'cell_mapping.json'
        with open(mapping_path, 'w') as f:
            json.dump(cell_mapping, f, indent=2)

        # Add to report
        self.registration_report['samples'][sample_id] = {
            'total_cells': len(cells),
            'cells_mapped': sum(1 for c in cell_mapping['cells'].values() if c['status'] == 'ALIGNED'),
            'mapping_path': str(mapping_path)
        }

    def _save_report(self):
        """Save registration report"""
        report_path = self.output_dir / 'registration_report.json'

        self.registration_report['summary'] = {
            'samples_processed': self.stats['samples_processed'],
            'samples_failed': self.stats['samples_failed'],
            'total_cells': self.stats['total_cells'],
            'cells_mapped': self.stats['cells_mapped'],
            'cells_with_issues': self.stats['cells_with_issues'],
            'mapping_success_rate': (self.stats['cells_mapped'] / self.stats['total_cells'] * 100)
                if self.stats['total_cells'] > 0 else 0
        }

        with open(report_path, 'w') as f:
            json.dump(self.registration_report, f, indent=2)

    def _print_results(self):
        """Print registration results"""
        summary = self.registration_report['summary']

        click.echo(f"\n✓ Stage 3 completed")
        click.echo(f"  → Output: {self.output_dir}")
        click.echo(f"\n  Samples processed: {summary['samples_processed']}")
        click.echo(f"  Samples failed: {summary['samples_failed']}")
        click.echo(f"  Total cells: {summary['total_cells']}")
        click.echo(f"  Cells mapped: {summary['cells_mapped']}")
        click.echo(f"  Mapping success: {summary['mapping_success_rate']:.1f}%")
        click.echo(f"  Cells with issues: {summary['cells_with_issues']}")

        if len(self.registration_report['errors']) > 0:
            click.echo(f"\n⚠ {len(self.registration_report['errors'])} errors. Check registration_report.json")

        if len(self.registration_report['warnings']) > 0:
            click.echo(f"⚠ {len(self.registration_report['warnings'])} warnings. Check registration_report.json")


def run_stage_3(config, output_dir):
    """Run Stage 3: Cell Registration/Mapping"""
    registrar = CellRegistrar(config, output_dir)
    return registrar.run()
