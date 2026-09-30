"""
Stage 1: Data Loading & Validation

Purpose:
    Load images, validate structure, create manifest of all files.

Input:
    - Raw image files in test_data/
    - Configuration file

Output:
    - outputs/01_validated/manifest.json - inventory of all images
    - outputs/01_validated/validation_report.json - file-level issues
"""

import json
import cv2
import numpy as np
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
import click


class DataValidator:
    """Validates microscopy image data"""

    def __init__(self, config, experiment_meta, output_dir):
        self.config = config
        self.experiment_meta = experiment_meta
        self.output_dir = Path(output_dir)
        self.manifest = {
            "created_at": datetime.now().isoformat(),
            "samples": {},
            "summary": {
                "total_samples": 0,
                "total_files": 0,
                "valid_samples": 0,
                "invalid_samples": 0
            }
        }
        self.validation_report = {
            "created_at": datetime.now().isoformat(),
            "errors": [],
            "warnings": [],
            "info": []
        }

    def run(self):
        """Execute Stage 1: Data Loading & Validation"""
        click.echo("\n" + "=" * 60)
        click.echo("[Stage 1] Data Loading & Validation")
        click.echo("=" * 60)

        # Load data paths
        data_dir = Path(self.experiment_meta['data_directory'])
        channels = self.experiment_meta['channels']

        # Scan for all samples
        click.echo("\n→ Scanning for image files...")
        samples = self._scan_samples(data_dir, channels)

        # Validate each sample
        click.echo(f"\n→ Validating {len(samples)} samples...")
        for sample_id, sample_data in tqdm(samples.items(), desc="Validating"):
            self._validate_sample(sample_id, sample_data)

        # Generate summary
        self._generate_summary()

        # Save outputs
        self._save_outputs()

        # Print results
        self._print_results()

        return self.manifest, self.validation_report

    def _scan_samples(self, data_dir, channels):
        """Scan data directory for all samples

        Supports two directory structures:
        1. Legacy: data_dir/Controls/Single channels/, data_dir/Treatment/Single channels/
        2. Experiment 1: data_dir/{condition}/ (no subdirectories)
        """
        samples = {}

        # Scan all subdirectories in data_dir (each is a condition/group)
        if not data_dir.exists():
            click.echo(f"⚠ Warning: Data directory not found: {data_dir}")
            return samples

        for condition_dir in data_dir.iterdir():
            if not condition_dir.is_dir():
                continue

            # Skip hidden directories
            if condition_dir.name.startswith('.'):
                continue

            # Determine group name (sanitize for use as directory name)
            group = condition_dir.name

            # Check if this is legacy structure (has "Single channels" subdirectory)
            single_channels_dir = condition_dir / "Single channels"
            if single_channels_dir.exists():
                # Legacy structure
                scan_dir = single_channels_dir
            else:
                # Experiment 1 structure (files directly in condition directory)
                scan_dir = condition_dir

            # Scan for channel files
            for channel in channels:
                pattern = f"{channel}-*.tif"
                for img_path in scan_dir.glob(pattern):
                    # Extract sample name from filename
                    # Format: C1-Untreated_1.tif → sample_name = "Untreated_1"
                    # Format: C2-4h 2-BP_5.tif → sample_name = "4h 2-BP_5"
                    sample_name = img_path.stem.split('-', 1)[1] if '-' in img_path.stem else img_path.stem

                    # Create unique sample_id combining group and sample_name
                    sample_id = f"{group}_{sample_name}"

                    if sample_id not in samples:
                        samples[sample_id] = {
                            "sample_id": sample_id,
                            "group": group,
                            "sample_name": sample_name,
                            "channels": {}
                        }

                    samples[sample_id]["channels"][channel] = {
                        "path": str(img_path),
                        "exists": True
                    }

        return samples

    def _validate_sample(self, sample_id, sample_data):
        """Validate a single sample"""
        channels_expected = set(self.experiment_meta['channels'])
        channels_found = set(sample_data['channels'].keys())

        # Check for missing channels
        missing_channels = channels_expected - channels_found
        if missing_channels:
            self.validation_report['errors'].append({
                "code": "E102",
                "level": "ERROR",
                "sample": sample_id,
                "message": f"Missing channels: {', '.join(missing_channels)}"
            })
            sample_data['valid'] = False
            return

        # Validate each channel image
        dimensions = {}
        for channel, channel_data in sample_data['channels'].items():
            img_path = Path(channel_data['path'])

            # Check file exists
            if not img_path.exists():
                self.validation_report['errors'].append({
                    "code": "E101",
                    "level": "ERROR",
                    "sample": sample_id,
                    "channel": channel,
                    "message": f"File not found: {img_path}"
                })
                sample_data['valid'] = False
                return

            # Try to load image
            try:
                img = cv2.imread(str(img_path), cv2.IMREAD_UNCHANGED)
                if img is None:
                    raise ValueError("Failed to load image")

                # Store image properties
                channel_data['shape'] = img.shape
                channel_data['dtype'] = str(img.dtype)
                channel_data['size_bytes'] = img_path.stat().st_size

                dimensions[channel] = img.shape

                # Check if image is blank
                if np.max(img) == 0:
                    self.validation_report['warnings'].append({
                        "code": "E201",
                        "level": "WARNING",
                        "sample": sample_id,
                        "channel": channel,
                        "message": "Image is completely black"
                    })

                # Check for saturation
                if img.dtype == np.uint8:
                    max_val = 255
                elif img.dtype == np.uint16:
                    max_val = 65535
                else:
                    max_val = np.max(img)

                saturated_pixels = np.sum(img == max_val)
                total_pixels = img.size
                saturation_percent = (saturated_pixels / total_pixels) * 100

                if saturation_percent > 30:
                    self.validation_report['warnings'].append({
                        "code": "E202",
                        "level": "WARNING",
                        "sample": sample_id,
                        "channel": channel,
                        "message": f"High saturation: {saturation_percent:.1f}% of pixels"
                    })

                # Check contrast
                img_range = np.max(img) - np.min(img)
                if img_range < 10:
                    self.validation_report['warnings'].append({
                        "code": "E203",
                        "level": "WARNING",
                        "sample": sample_id,
                        "channel": channel,
                        "message": f"Low contrast: range = {img_range}"
                    })

            except Exception as e:
                self.validation_report['errors'].append({
                    "code": "E103",
                    "level": "ERROR",
                    "sample": sample_id,
                    "channel": channel,
                    "message": f"Corrupted file: {str(e)}"
                })
                sample_data['valid'] = False
                return

        # Check dimension consistency across channels
        shapes = list(dimensions.values())
        if len(set(shapes)) > 1:
            self.validation_report['errors'].append({
                "code": "E104",
                "level": "ERROR",
                "sample": sample_id,
                "message": f"Dimension mismatch across channels: {dimensions}"
            })
            sample_data['valid'] = False
            return

        # Sample is valid
        sample_data['valid'] = True
        self.manifest['samples'][sample_id] = sample_data

    def _generate_summary(self):
        """Generate validation summary"""
        total_samples = len(self.manifest['samples'])
        valid_samples = sum(1 for s in self.manifest['samples'].values() if s.get('valid', False))
        invalid_samples = total_samples - valid_samples

        total_files = sum(
            len(s['channels'])
            for s in self.manifest['samples'].values()
        )

        self.manifest['summary'] = {
            "total_samples": total_samples,
            "total_files": total_files,
            "valid_samples": valid_samples,
            "invalid_samples": invalid_samples,
            "error_count": len(self.validation_report['errors']),
            "warning_count": len(self.validation_report['warnings'])
        }

        self.validation_report['summary'] = {
            "total_errors": len(self.validation_report['errors']),
            "total_warnings": len(self.validation_report['warnings']),
            "total_info": len(self.validation_report['info']),
            "samples_with_errors": len(set(
                e['sample'] for e in self.validation_report['errors'] if 'sample' in e
            )),
            "samples_with_warnings": len(set(
                w['sample'] for w in self.validation_report['warnings'] if 'sample' in w
            ))
        }

    def _save_outputs(self):
        """Save manifest and validation report"""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Save manifest
        manifest_path = self.output_dir / "manifest.json"
        with open(manifest_path, 'w') as f:
            json.dump(self.manifest, f, indent=2)

        # Save validation report
        report_path = self.output_dir / "validation_report.json"
        with open(report_path, 'w') as f:
            json.dump(self.validation_report, f, indent=2)

    def _print_results(self):
        """Print validation results"""
        summary = self.manifest['summary']

        click.echo(f"\n✓ Stage 1 completed")
        click.echo(f"  → Output: {self.output_dir}")
        click.echo(f"\n  Samples: {summary['valid_samples']}/{summary['total_samples']} valid")
        click.echo(f"  Files: {summary['total_files']} total")
        click.echo(f"  Errors: {summary['error_count']}")
        click.echo(f"  Warnings: {summary['warning_count']}")

        if summary['error_count'] > 0:
            click.echo(f"\n⚠ {summary['error_count']} errors found. Check validation_report.json")

        if summary['warning_count'] > 0:
            click.echo(f"⚠ {summary['warning_count']} warnings found. Check validation_report.json")


def run_stage_1(config, experiment_meta, output_dir):
    """Run Stage 1: Data Loading & Validation"""
    validator = DataValidator(config, experiment_meta, output_dir)
    return validator.run()
