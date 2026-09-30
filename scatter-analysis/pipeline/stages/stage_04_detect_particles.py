"""
Stage 4: Particle Detection

Purpose:
    Identify individual particles within each cell for all channels.

Input:
    - outputs/02_segmented/[group]/[sample_id]/[channel]/cell_[n].tif
    - outputs/03_registered/[group]/[sample_id]/cell_mapping.json

Output:
    - outputs/04_particles/[group]/[sample_id]/[channel]/cell_[n]_particles.json

Example output:
{
  "cell_id": "Control_1_cell_1",
  "channel": "C1",
  "total_particles": 145,
  "cell_center": [150, 190],
  "particles": [
    {
      "id": 0,
      "position": [145, 185],
      "area": 12,
      "intensity_mean": 187.3,
      "intensity_max": 255,
      "distance_from_cell_center": 5.1
    },
    ...
  ]
}
"""

import json
import cv2
import numpy as np
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
import click


class ParticleDetector:
    """Detects particles within cells"""

    def __init__(self, config, output_dir):
        self.config = config
        self.output_dir = Path(output_dir)
        self.segmented_dir = Path("outputs/02_segmented")
        self.registered_dir = Path("outputs/03_registered")

        # Particle detection parameters from config
        particle_config = config.get('particle_detection', {})
        self.threshold_method = particle_config.get('threshold_method', 'otsu')
        self.min_particle_size = particle_config.get('min_particle_size', 1)
        self.max_particle_size = particle_config.get('max_particle_size', 1000)
        self.blur_sigma = particle_config.get('blur_sigma', 1.0)

        # Statistics
        self.stats = {
            'samples_processed': 0,
            'cells_processed': 0,
            'total_particles': 0,
            'particles_by_channel': {'C1': 0, 'C2': 0, 'C3': 0}
        }

        self.detection_report = {
            'created_at': datetime.now().isoformat(),
            'samples': {},
            'errors': [],
            'warnings': []
        }

    def run(self):
        """Execute Stage 4: Particle Detection"""
        click.echo("\n" + "=" * 60)
        click.echo("[Stage 4] Particle Detection")
        click.echo("=" * 60)

        # Load registration report to get list of samples
        reg_report_path = self.registered_dir / "registration_report.json"
        with open(reg_report_path, 'r') as f:
            reg_report = json.load(f)

        samples = reg_report.get('samples', {})

        click.echo(f"\n→ Detecting particles in cells for {len(samples)} samples...")

        for sample_id, sample_data in tqdm(samples.items(), desc="Detecting particles"):
            try:
                self._process_sample(sample_id)
                self.stats['samples_processed'] += 1
            except Exception as e:
                click.echo(f"\n✗ Error processing {sample_id}: {str(e)}")
                self.detection_report['errors'].append({
                    'sample': sample_id,
                    'error': str(e)
                })

        # Save detection report
        self._save_report()

        # Print results
        self._print_results()

        return self.detection_report

    def _process_sample(self, sample_id):
        """Process all cells in a sample"""

        # Parse sample_id to get group and sample_name
        parts = sample_id.split('_', 1)
        group = parts[0]
        sample_name = parts[1]

        # Load cell mapping
        mapping_file = self.registered_dir / group / sample_name / 'cell_mapping.json'
        with open(mapping_file, 'r') as f:
            mapping = json.load(f)

        # Create output directories
        sample_output_dir = self.output_dir / group / sample_name
        for channel in ['C1', 'C2', 'C3']:
            (sample_output_dir / channel).mkdir(parents=True, exist_ok=True)

        sample_particle_counts = {
            'sample_id': sample_id,
            'cells': {},
            'total_particles': 0
        }

        # Process each cell
        for cell_id, cell_info in mapping['cells'].items():
            if cell_info['status'] != 'ALIGNED':
                continue

            # Load cell metadata
            metadata_file = Path(cell_info['metadata_path'])
            with open(metadata_file, 'r') as f:
                cell_metadata = json.load(f)

            # Detect particles in each channel
            for channel, ch_data in cell_metadata['channels'].items():
                particles = self._detect_particles_in_cell(
                    cell_id, channel, ch_data, sample_id
                )

                # Save particle data
                output_file = sample_output_dir / channel / f"{cell_id}_particles.json"
                with open(output_file, 'w') as f:
                    json.dump(particles, f, indent=2)

                # Update statistics
                particle_count = particles['total_particles']
                self.stats['total_particles'] += particle_count
                self.stats['particles_by_channel'][channel] += particle_count
                sample_particle_counts['total_particles'] += particle_count

                if cell_id not in sample_particle_counts['cells']:
                    sample_particle_counts['cells'][cell_id] = {}
                sample_particle_counts['cells'][cell_id][channel] = particle_count

            self.stats['cells_processed'] += 1

        self.detection_report['samples'][sample_id] = sample_particle_counts

    def _detect_particles_in_cell(self, cell_id, channel, ch_data, sample_id):
        """Detect particles in a single cell image"""

        # Load cell image
        cell_img_path = Path(ch_data['path'])
        img = cv2.imread(str(cell_img_path), cv2.IMREAD_UNCHANGED)

        if img is None:
            raise ValueError(f"Failed to load cell image: {cell_img_path}")

        # Convert to grayscale if RGB
        if len(img.shape) == 3:
            gray = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
        else:
            gray = img.copy()

        # Preprocessing
        # 1. Gaussian blur to reduce noise
        kernel_size = int(2 * round(2 * self.blur_sigma) + 1)  # Ensure odd kernel
        blurred = cv2.GaussianBlur(gray, (kernel_size, kernel_size), self.blur_sigma)

        # 2. Thresholding
        if self.threshold_method == 'otsu':
            threshold_value, binary = cv2.threshold(
                blurred, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU
            )
        else:
            # Default to Otsu
            threshold_value, binary = cv2.threshold(
                blurred, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU
            )

        # 3. Find connected components (particles)
        num_labels, labels, stats, centroids = cv2.connectedComponentsWithStats(
            binary, connectivity=8
        )

        # Calculate cell center (center of the image for now)
        cell_height, cell_width = gray.shape
        cell_center = [cell_width // 2, cell_height // 2]

        # Extract particle information
        particles = []

        # Skip label 0 (background)
        for i in range(1, num_labels):
            area = stats[i, cv2.CC_STAT_AREA]

            # Filter by size
            if area < self.min_particle_size or area > self.max_particle_size:
                continue

            # Get particle position (centroid)
            cx, cy = centroids[i]
            position = [float(cx), float(cy)]

            # Calculate distance from cell center
            distance = np.sqrt((cx - cell_center[0])**2 + (cy - cell_center[1])**2)

            # Get intensity statistics for this particle
            mask = (labels == i).astype(np.uint8)
            particle_pixels = gray[mask == 1]

            if len(particle_pixels) > 0:
                intensity_mean = float(np.mean(particle_pixels))
                intensity_max = float(np.max(particle_pixels))
                intensity_min = float(np.min(particle_pixels))
                intensity_std = float(np.std(particle_pixels))
            else:
                intensity_mean = 0.0
                intensity_max = 0.0
                intensity_min = 0.0
                intensity_std = 0.0

            particles.append({
                'id': i - 1,  # 0-indexed
                'position': position,
                'area': int(area),
                'distance_from_cell_center': float(distance),
                'intensity_mean': intensity_mean,
                'intensity_max': intensity_max,
                'intensity_min': intensity_min,
                'intensity_std': intensity_std
            })

        # Check for potential issues
        if len(particles) == 0:
            self.detection_report['warnings'].append({
                'code': 'E501',
                'sample': sample_id,
                'cell': cell_id,
                'channel': channel,
                'message': 'No particles detected'
            })
        elif len(particles) > 50000:
            self.detection_report['warnings'].append({
                'code': 'E502',
                'sample': sample_id,
                'cell': cell_id,
                'channel': channel,
                'message': f'Very high particle count: {len(particles)}'
            })

        return {
            'cell_id': f"{sample_id}_{cell_id}",
            'channel': channel,
            'total_particles': len(particles),
            'cell_center': cell_center,
            'cell_shape': [cell_height, cell_width],
            'threshold_value': float(threshold_value),
            'particles': particles
        }

    def _save_report(self):
        """Save detection report"""
        report_path = self.output_dir / 'particle_detection_report.json'

        avg_particles_per_cell = (self.stats['total_particles'] / self.stats['cells_processed']
                                   if self.stats['cells_processed'] > 0 else 0)

        self.detection_report['summary'] = {
            'samples_processed': self.stats['samples_processed'],
            'cells_processed': self.stats['cells_processed'],
            'total_particles': self.stats['total_particles'],
            'particles_by_channel': self.stats['particles_by_channel'],
            'avg_particles_per_cell': avg_particles_per_cell
        }

        with open(report_path, 'w') as f:
            json.dump(self.detection_report, f, indent=2)

    def _print_results(self):
        """Print detection results"""
        summary = self.detection_report['summary']

        click.echo(f"\n✓ Stage 4 completed")
        click.echo(f"  → Output: {self.output_dir}")
        click.echo(f"\n  Samples processed: {summary['samples_processed']}")
        click.echo(f"  Cells processed: {summary['cells_processed']}")
        click.echo(f"  Total particles: {summary['total_particles']:,}")
        click.echo(f"  Avg particles/cell: {summary['avg_particles_per_cell']:.1f}")
        click.echo(f"\n  Particles by channel:")
        for ch, count in summary['particles_by_channel'].items():
            click.echo(f"    {ch}: {count:,}")

        if len(self.detection_report['warnings']) > 0:
            click.echo(f"\n⚠ {len(self.detection_report['warnings'])} warnings. Check particle_detection_report.json")


def run_stage_4(config, output_dir):
    """Run Stage 4: Particle Detection"""
    detector = ParticleDetector(config, output_dir)
    return detector.run()
