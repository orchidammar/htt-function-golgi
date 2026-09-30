"""
Stage 5: Scatter Metrics Calculation

Purpose:
    Calculate quantitative metrics describing particle scattering/spread.

Primary Metric: MEDIAN distance from center (outlier-resistant)

Input:
    - outputs/04_particles/[group]/[sample_id]/[channel]/cell_[n]_particles.json

Output:
    - outputs/05_metrics/[group]/[sample_id]/[channel]/cell_[n]_metrics.json

Metrics Categories:
1. Distance-Based (Primary): mean, median, std, max, min, range
2. Area-Based: convex_hull_area, occupied_area, percent_cell_occupied
3. Density-Based: particle_density, particles_per_unit_area
4. Dispersion: coefficient_of_variation, quartile_dispersion
5. Percentiles: p10, p25, p50, p75, p90, p95, p99
"""

import json
import numpy as np
from pathlib import Path
from datetime import datetime
from scipy.spatial import ConvexHull
from scipy.stats import pearsonr, spearmanr
from tqdm import tqdm
import click
import cv2


class MetricsCalculator:
    """Calculates scatter metrics for particles"""

    def __init__(self, config, output_dir):
        self.config = config
        self.output_dir = Path(output_dir)
        self.particles_dir = Path("outputs/04_particles")

        # Metrics config
        metrics_config = config.get('metrics', {})
        self.primary_metric = metrics_config.get('primary_metric', 'median_distance')

        # Statistics
        self.stats = {
            'samples_processed': 0,
            'cells_processed': 0,
            'cells_skipped': 0
        }

        self.metrics_report = {
            'created_at': datetime.now().isoformat(),
            'samples': {},
            'errors': [],
            'warnings': []
        }

    def run(self):
        """Execute Stage 5: Scatter Metrics Calculation"""
        click.echo("\n" + "=" * 60)
        click.echo("[Stage 5] Scatter Metrics Calculation")
        click.echo("=" * 60)

        # Load particle detection report
        particle_report_path = self.particles_dir / "particle_detection_report.json"
        with open(particle_report_path, 'r') as f:
            particle_report = json.load(f)

        samples = particle_report.get('samples', {})

        click.echo(f"\n→ Calculating scatter metrics for {len(samples)} samples...")
        click.echo(f"   Primary metric: {self.primary_metric}")

        for sample_id, sample_data in tqdm(samples.items(), desc="Calculating metrics"):
            try:
                self._process_sample(sample_id, sample_data)
                self.stats['samples_processed'] += 1
            except Exception as e:
                click.echo(f"\n✗ Error processing {sample_id}: {str(e)}")
                self.metrics_report['errors'].append({
                    'sample': sample_id,
                    'error': str(e)
                })

        # Save metrics report
        self._save_report()

        # Print results
        self._print_results()

        return self.metrics_report

    def _process_sample(self, sample_id, sample_data):
        """Process all cells in a sample"""

        try:
            # Parse sample_id
            parts = sample_id.split('_', 1)
            group = parts[0]
            sample_name = parts[1]
        except Exception as e:
            raise Exception(f"Error parsing sample_id '{sample_id}': {e}")

        # Create output directories
        sample_output_dir = self.output_dir / group / sample_name
        for channel in ['C1', 'C2', 'C3']:
            (sample_output_dir / channel).mkdir(parents=True, exist_ok=True)

        sample_metrics_summary = {
            'sample_id': sample_id,
            'cells': {}
        }

        # Process each cell
        for cell_id, cell_particle_counts in sample_data.get('cells', {}).items():
            for channel in ['C1', 'C2', 'C3']:
                particle_count = cell_particle_counts.get(channel, 0)

                # Load particle data
                particle_file = self.particles_dir / group / sample_name / channel / f"{cell_id}_particles.json"

                if not particle_file.exists():
                    continue

                with open(particle_file, 'r') as f:
                    particle_data = json.load(f)

                # Calculate metrics
                metrics = self._calculate_metrics(particle_data, sample_id, cell_id, channel)

                # Save metrics
                output_file = sample_output_dir / channel / f"{cell_id}_metrics.json"
                with open(output_file, 'w') as f:
                    json.dump(metrics, f, indent=2)

                # Add to summary
                if cell_id not in sample_metrics_summary['cells']:
                    sample_metrics_summary['cells'][cell_id] = {}
                sample_metrics_summary['cells'][cell_id][channel] = {
                    'median_distance': metrics.get('median_distance_from_center'),
                    'mean_distance': metrics.get('mean_distance_from_center'),
                    'particle_count': metrics.get('total_particles')
                }

                self.stats['cells_processed'] += 1

            # Calculate intensity-based cross-channel correlation (C1-C2 co-scattering)
            # This uses pixel-by-pixel analysis instead of particle detection
            correlation_metrics = self._calculate_intensity_based_correlation(group, sample_name, cell_id)

            if correlation_metrics and correlation_metrics.get('status') == 'SUCCESS':
                # Save correlation metrics to separate file
                correlation_output_dir = sample_output_dir / 'correlation'
                correlation_output_dir.mkdir(parents=True, exist_ok=True)
                correlation_file = correlation_output_dir / f"{cell_id}_correlation.json"

                with open(correlation_file, 'w') as f:
                    json.dump(correlation_metrics, f, indent=2)

                # Add to summary
                if cell_id in sample_metrics_summary['cells']:
                    sample_metrics_summary['cells'][cell_id]['correlation'] = {
                        'pixel_pearson_r': correlation_metrics.get('pixel_pearson_r'),
                        'manders_m1': correlation_metrics.get('manders_m1'),
                        'manders_m2': correlation_metrics.get('manders_m2'),
                        'intensity_centroid_distance': correlation_metrics.get('intensity_centroid_distance'),
                        'c1_scatter_index': correlation_metrics.get('c1_scatter_index'),
                        'c2_scatter_index': correlation_metrics.get('c2_scatter_index'),
                        'interpretation': correlation_metrics.get('interpretation')
                    }

        self.metrics_report['samples'][sample_id] = sample_metrics_summary

    def _calculate_metrics(self, particle_data, sample_id, cell_id, channel):
        """Calculate all scatter metrics for a cell"""

        particles = particle_data.get('particles', [])
        total_particles = len(particles)
        cell_center = particle_data.get('cell_center', [0, 0])
        cell_shape = particle_data.get('cell_shape', [0, 0])
        cell_area = cell_shape[0] * cell_shape[1] if len(cell_shape) >= 2 else 0

        metrics = {
            'cell_id': particle_data.get('cell_id'),
            'channel': channel,
            'total_particles': total_particles,
            'cell_center': cell_center,
            'cell_area': cell_area,
            'threshold_value': particle_data.get('threshold_value')
        }

        # Handle case with no particles
        if total_particles == 0:
            self.stats['cells_skipped'] += 1
            metrics['status'] = 'NO_PARTICLES'
            return metrics

        # Extract particle positions and distances
        positions = np.array([p['position'] for p in particles])
        distances = np.array([p['distance_from_cell_center'] for p in particles])
        areas = np.array([p['area'] for p in particles])
        intensities = np.array([p['intensity_mean'] for p in particles])

        # 1. DISTANCE-BASED METRICS (Primary)
        metrics['mean_distance_from_center'] = float(np.mean(distances))
        metrics['median_distance_from_center'] = float(np.median(distances))  # PRIMARY METRIC
        metrics['std_distance_from_center'] = float(np.std(distances))
        metrics['max_distance_from_center'] = float(np.max(distances))
        metrics['min_distance_from_center'] = float(np.min(distances))
        metrics['distance_range'] = float(np.max(distances) - np.min(distances))

        # Percentiles
        percentiles = [10, 25, 50, 75, 90, 95, 99]
        for p in percentiles:
            if total_particles > 0:
                metrics[f'distance_p{p}'] = float(np.percentile(distances, p))
            else:
                metrics[f'distance_p{p}'] = 0.0

        # 2. AREA-BASED METRICS
        # Convex hull area (spatial extent)
        if total_particles >= 3:
            try:
                hull = ConvexHull(positions)
                metrics['convex_hull_area'] = float(hull.volume)  # 2D area
                metrics['convex_hull_to_cell_ratio'] = float(hull.volume / cell_area) if cell_area > 0 else 0.0
            except Exception as e:
                metrics['convex_hull_area'] = 0.0
                metrics['convex_hull_to_cell_ratio'] = 0.0
        else:
            metrics['convex_hull_area'] = 0.0
            metrics['convex_hull_to_cell_ratio'] = 0.0

        # Occupied area (sum of particle areas)
        metrics['occupied_area'] = float(np.sum(areas))
        metrics['percent_cell_occupied'] = float((np.sum(areas) / cell_area) * 100) if cell_area > 0 else 0.0

        # 3. DENSITY-BASED METRICS
        metrics['particle_density'] = float(total_particles / cell_area) if cell_area > 0 else 0.0

        if metrics['convex_hull_area'] > 0:
            metrics['particle_density_in_hull'] = float(total_particles / metrics['convex_hull_area'])
        else:
            metrics['particle_density_in_hull'] = 0.0

        # 4. DISPERSION METRICS
        if metrics['mean_distance_from_center'] > 0:
            metrics['coefficient_of_variation'] = float(metrics['std_distance_from_center'] / metrics['mean_distance_from_center'])
        else:
            metrics['coefficient_of_variation'] = 0.0

        # IQR-based dispersion
        if total_particles > 1:
            q75, q25 = np.percentile(distances, [75, 25])
            iqr = q75 - q25
            metrics['quartile_dispersion'] = float(iqr)
        else:
            metrics['quartile_dispersion'] = 0.0

        # 5. OUTLIER DETECTION (IQR method)
        if total_particles > 3:
            q75, q25 = np.percentile(distances, [75, 25])
            iqr = q75 - q25
            lower_bound = q25 - 1.5 * iqr
            upper_bound = q75 + 1.5 * iqr
            outliers = np.sum((distances < lower_bound) | (distances > upper_bound))
            outlier_percent = (outliers / total_particles) * 100

            metrics['outliers_detected'] = int(outliers)
            metrics['outlier_percent'] = float(outlier_percent)
            metrics['metric_used'] = 'MEDIAN' if outlier_percent > 5 else 'MEAN'
        else:
            metrics['outliers_detected'] = 0
            metrics['outlier_percent'] = 0.0
            metrics['metric_used'] = 'MEAN'

        # 6. INTENSITY METRICS
        metrics['mean_particle_intensity'] = float(np.mean(intensities))
        metrics['median_particle_intensity'] = float(np.median(intensities))
        metrics['total_intensity'] = float(np.sum(intensities))

        # 7. SHAPE METRICS (aspect ratio of distribution)
        if total_particles >= 2:
            x_coords = positions[:, 0]
            y_coords = positions[:, 1]
            x_range = np.max(x_coords) - np.min(x_coords)
            y_range = np.max(y_coords) - np.min(y_coords)

            if min(x_range, y_range) > 0:
                metrics['distribution_aspect_ratio'] = float(max(x_range, y_range) / min(x_range, y_range))
            else:
                metrics['distribution_aspect_ratio'] = 1.0
        else:
            metrics['distribution_aspect_ratio'] = 1.0

        metrics['status'] = 'SUCCESS'

        return metrics

    def _calculate_cross_channel_correlation(self, group, sample_name, cell_id):
        """Calculate C1-C2 spatial correlation metrics for co-scattering analysis

        This measures whether protein X (C2/GREEN) and structure Y (C1/RED) scatter together.

        NOTE: C3 (BLUE) is the central structure used for positioning only, NOT analyzed.

        Returns correlation metrics or None if data unavailable.
        """

        # Load C1 (RED/Structure Y) and C2 (GREEN/Protein X) particle data
        c1_file = self.particles_dir / group / sample_name / 'C1' / f"{cell_id}_particles.json"
        c2_file = self.particles_dir / group / sample_name / 'C2' / f"{cell_id}_particles.json"

        if not (c1_file.exists() and c2_file.exists()):
            return {'status': 'MISSING_DATA', 'c1_exists': c1_file.exists(), 'c2_exists': c2_file.exists()}

        # Load particle data
        with open(c1_file, 'r') as f:
            c1_data = json.load(f)
        with open(c2_file, 'r') as f:
            c2_data = json.load(f)

        c1_particles = c1_data.get('particles', [])  # Structure Y (RED)
        c2_particles = c2_data.get('particles', [])  # Protein X (GREEN)

        # Initialize metrics
        metrics = {
            'cell_id': cell_id,
            'c1_particle_count': len(c1_particles),  # Structure Y
            'c2_particle_count': len(c2_particles)   # Protein X
        }

        # Check for particles
        if len(c1_particles) == 0 or len(c2_particles) == 0:
            metrics['status'] = 'NO_PARTICLES'
            return metrics

        # Need at least 3 particles in each channel for meaningful correlation
        if len(c1_particles) < 3 or len(c2_particles) < 3:
            metrics['status'] = 'INSUFFICIENT_PARTICLES'
            return metrics

        # Extract distances from cell center
        c1_distances = np.array([p['distance_from_cell_center'] for p in c1_particles])  # Structure Y
        c2_distances = np.array([p['distance_from_cell_center'] for p in c2_particles])  # Protein X

        # Extract positions
        c1_positions = np.array([p['position'] for p in c1_particles])  # Structure Y
        c2_positions = np.array([p['position'] for p in c2_particles])  # Protein X

        # 1. PEARSON CORRELATION (Linear correlation of scatter distances)
        # Measures if C1 (Structure Y) and C2 (Protein X) have similar scatter magnitude distributions
        try:
            # Correlation between distance distributions
            # Create distance bins and correlate histograms
            max_dist = max(c1_distances.max(), c2_distances.max())
            bins = np.linspace(0, max_dist, 20)
            c1_hist, _ = np.histogram(c1_distances, bins=bins)
            c2_hist, _ = np.histogram(c2_distances, bins=bins)

            if c1_hist.sum() > 0 and c2_hist.sum() > 0:
                r, p_val = pearsonr(c1_hist, c2_hist)
                metrics['pearson_r_distance_distribution'] = float(r)
                metrics['pearson_p_distance_distribution'] = float(p_val)
            else:
                metrics['pearson_r_distance_distribution'] = None
                metrics['pearson_p_distance_distribution'] = None
        except Exception as e:
            metrics['pearson_r_distance_distribution'] = None
            metrics['pearson_p_distance_distribution'] = None
            metrics['pearson_error'] = str(e)

        # 2. SPEARMAN CORRELATION (Rank-based, robust to outliers)
        try:
            rho, p_val = spearmanr(c1_hist, c2_hist)
            metrics['spearman_rho_distance_distribution'] = float(rho)
            metrics['spearman_p_distance_distribution'] = float(p_val)
        except Exception as e:
            metrics['spearman_rho_distance_distribution'] = None
            metrics['spearman_p_distance_distribution'] = None
            metrics['spearman_error'] = str(e)

        # 3. SCATTER METRICS CORRELATION
        # Correlate summary statistics (more robust for different particle counts)
        c1_median = float(np.median(c1_distances))  # Structure Y
        c2_median = float(np.median(c2_distances))  # Protein X
        c1_mean = float(np.mean(c1_distances))
        c2_mean = float(np.mean(c2_distances))

        metrics['c1_median_distance'] = c1_median
        metrics['c2_median_distance'] = c2_median
        metrics['c1_mean_distance'] = c1_mean
        metrics['c2_mean_distance'] = c2_mean

        # Difference in scatter (how different are median distances?)
        metrics['median_distance_difference'] = abs(c1_median - c2_median)
        metrics['mean_distance_difference'] = abs(c1_mean - c2_mean)

        # 4. SPATIAL PROXIMITY (Centroid distance)
        c1_centroid = np.mean(c1_positions, axis=0)  # Structure Y
        c2_centroid = np.mean(c2_positions, axis=0)  # Protein X
        centroid_distance = np.linalg.norm(c1_centroid - c2_centroid)

        metrics['c1_centroid'] = c1_centroid.tolist()
        metrics['c2_centroid'] = c2_centroid.tolist()
        metrics['centroid_distance'] = float(centroid_distance)

        # 5. MANDERS COLOCALIZATION COEFFICIENTS (M1 and M2)
        # M1 = fraction of C2 (Protein X) intensity overlapping with C1 (Structure Y)
        # M2 = fraction of C1 (Structure Y) intensity overlapping with C2 (Protein X)
        # Use spatial proximity (within threshold distance)

        colocalization_threshold = 5.0  # pixels (adjust based on microscope resolution)

        # For each particle, check if any particle from other channel is nearby
        c1_intensities = np.array([p.get('intensity_mean', 1.0) for p in c1_particles])
        c2_intensities = np.array([p.get('intensity_mean', 1.0) for p in c2_particles])

        c1_coloc_intensity = 0.0
        c2_coloc_intensity = 0.0

        for i, c2_pos in enumerate(c2_positions):
            # Check if this C2 (Protein X) particle has a C1 (Structure Y) neighbor
            distances_to_c1 = np.linalg.norm(c1_positions - c2_pos, axis=1)
            if np.any(distances_to_c1 < colocalization_threshold):
                c2_coloc_intensity += c2_intensities[i]

        for j, c1_pos in enumerate(c1_positions):
            # Check if this C1 (Structure Y) particle has a C2 (Protein X) neighbor
            distances_to_c2 = np.linalg.norm(c2_positions - c1_pos, axis=1)
            if np.any(distances_to_c2 < colocalization_threshold):
                c1_coloc_intensity += c1_intensities[j]

        # Calculate Manders coefficients
        c1_total_intensity = np.sum(c1_intensities)
        c2_total_intensity = np.sum(c2_intensities)

        if c2_total_intensity > 0:
            metrics['manders_m1'] = float(c2_coloc_intensity / c2_total_intensity)
        else:
            metrics['manders_m1'] = 0.0

        if c1_total_intensity > 0:
            metrics['manders_m2'] = float(c1_coloc_intensity / c1_total_intensity)
        else:
            metrics['manders_m2'] = 0.0

        metrics['colocalization_threshold_pixels'] = colocalization_threshold

        # 6. OVERLAP INTERPRETATION
        # Interpret results based on standard thresholds
        pearson_r = metrics.get('pearson_r_distance_distribution')
        manders_avg = (metrics['manders_m1'] + metrics['manders_m2']) / 2.0

        if pearson_r is not None and pearson_r > 0.6 and manders_avg > 0.7:
            metrics['interpretation'] = 'HIGH_COSCATTERING'
        elif pearson_r is not None and pearson_r > 0.3 and manders_avg > 0.5:
            metrics['interpretation'] = 'MODERATE_COSCATTERING'
        elif pearson_r is not None and pearson_r < 0.3:
            metrics['interpretation'] = 'LOW_COSCATTERING'
        else:
            metrics['interpretation'] = 'INDETERMINATE'

        metrics['status'] = 'SUCCESS'

        return metrics

    def _calculate_intensity_based_correlation(self, group, sample_name, cell_id):
        """Calculate intensity-based (pixel-by-pixel) correlation between C1 and C2

        This analyzes the raw fluorescence intensities without particle detection,
        measuring:
        1. Pixel-wise correlation between C1 (RED/Structure Y) and C2 (GREEN/Protein X)
        2. Intensity dispersion (how scattered bright pixels are)
        3. Radial intensity profiles (intensity vs distance from center)
        4. Intensity-based colocalization (Manders coefficients on full images)

        This works for ALL cells regardless of particle count.
        """

        # Get paths to segmented cell images
        segmented_dir = Path("outputs/02_segmented")
        c1_img_path = segmented_dir / group / sample_name / 'C1' / f"{cell_id}.tif"
        c2_img_path = segmented_dir / group / sample_name / 'C2' / f"{cell_id}.tif"

        if not (c1_img_path.exists() and c2_img_path.exists()):
            return {'status': 'MISSING_IMAGES',
                    'c1_exists': c1_img_path.exists(),
                    'c2_exists': c2_img_path.exists()}

        # Load images
        c1_img = cv2.imread(str(c1_img_path), cv2.IMREAD_UNCHANGED)
        c2_img = cv2.imread(str(c2_img_path), cv2.IMREAD_UNCHANGED)

        if c1_img is None or c2_img is None:
            return {'status': 'IMAGE_READ_ERROR'}

        # Convert to grayscale if needed
        if len(c1_img.shape) == 3:
            c1_img = cv2.cvtColor(c1_img, cv2.COLOR_BGR2GRAY)
        if len(c2_img.shape) == 3:
            c2_img = cv2.cvtColor(c2_img, cv2.COLOR_BGR2GRAY)

        # Ensure images are the same size
        if c1_img.shape != c2_img.shape:
            return {'status': 'SIZE_MISMATCH',
                    'c1_shape': c1_img.shape,
                    'c2_shape': c2_img.shape}

        # Initialize metrics
        metrics = {
            'cell_id': cell_id,
            'image_shape': list(c1_img.shape)
        }

        # Convert to float for calculations
        c1_float = c1_img.astype(np.float64)
        c2_float = c2_img.astype(np.float64)

        # Get cell center (use image center as approximation)
        height, width = c1_img.shape
        cell_center_y, cell_center_x = height / 2, width / 2
        metrics['cell_center'] = [cell_center_x, cell_center_y]

        # Create coordinate grids for distance calculations
        y_coords, x_coords = np.mgrid[0:height, 0:width]
        distances = np.sqrt((x_coords - cell_center_x)**2 + (y_coords - cell_center_y)**2)

        # 1. PIXEL-WISE CORRELATION
        # Flatten images to 1D arrays for correlation
        c1_flat = c1_float.flatten()
        c2_flat = c2_float.flatten()

        # Only correlate non-zero pixels (ignore background)
        mask = (c1_flat > 0) | (c2_flat > 0)

        if mask.sum() < 10:  # Need at least 10 non-zero pixels
            metrics['status'] = 'INSUFFICIENT_SIGNAL'
            return metrics

        c1_nonzero = c1_flat[mask]
        c2_nonzero = c2_flat[mask]

        try:
            pearson_r, pearson_p = pearsonr(c1_nonzero, c2_nonzero)
            metrics['pixel_pearson_r'] = float(pearson_r)
            metrics['pixel_pearson_p'] = float(pearson_p)
        except Exception as e:
            metrics['pixel_pearson_r'] = None
            metrics['pixel_pearson_p'] = None
            metrics['pearson_error'] = str(e)

        try:
            spearman_rho, spearman_p = spearmanr(c1_nonzero, c2_nonzero)
            metrics['pixel_spearman_rho'] = float(spearman_rho)
            metrics['pixel_spearman_p'] = float(spearman_p)
        except Exception as e:
            metrics['pixel_spearman_rho'] = None
            metrics['pixel_spearman_p'] = None
            metrics['spearman_error'] = str(e)

        # 2. INTENSITY STATISTICS
        metrics['c1_total_intensity'] = float(c1_float.sum())
        metrics['c2_total_intensity'] = float(c2_float.sum())
        metrics['c1_mean_intensity'] = float(c1_nonzero.mean()) if len(c1_nonzero) > 0 else 0.0
        metrics['c2_mean_intensity'] = float(c2_nonzero.mean()) if len(c2_nonzero) > 0 else 0.0
        metrics['c1_max_intensity'] = float(c1_float.max())
        metrics['c2_max_intensity'] = float(c2_float.max())
        metrics['c1_nonzero_pixel_count'] = int((c1_float > 0).sum())
        metrics['c2_nonzero_pixel_count'] = int((c2_float > 0).sum())

        # 3. INTENSITY-WEIGHTED CENTROIDS (center of mass)
        if metrics['c1_total_intensity'] > 0:
            c1_centroid_x = float((c1_float * x_coords).sum() / metrics['c1_total_intensity'])
            c1_centroid_y = float((c1_float * y_coords).sum() / metrics['c1_total_intensity'])
            metrics['c1_intensity_centroid'] = [c1_centroid_x, c1_centroid_y]
        else:
            metrics['c1_intensity_centroid'] = [cell_center_x, cell_center_y]

        if metrics['c2_total_intensity'] > 0:
            c2_centroid_x = float((c2_float * x_coords).sum() / metrics['c2_total_intensity'])
            c2_centroid_y = float((c2_float * y_coords).sum() / metrics['c2_total_intensity'])
            metrics['c2_intensity_centroid'] = [c2_centroid_x, c2_centroid_y]
        else:
            metrics['c2_intensity_centroid'] = [cell_center_x, cell_center_y]

        # Distance between intensity centroids
        c1_cent = np.array(metrics['c1_intensity_centroid'])
        c2_cent = np.array(metrics['c2_intensity_centroid'])
        metrics['intensity_centroid_distance'] = float(np.linalg.norm(c1_cent - c2_cent))

        # 4. INTENSITY DISPERSION (scatter index)
        # Intensity-weighted standard deviation from cell center
        if metrics['c1_total_intensity'] > 0:
            c1_weighted_dist_sq = (c1_float * distances**2).sum() / metrics['c1_total_intensity']
            metrics['c1_intensity_weighted_std'] = float(np.sqrt(c1_weighted_dist_sq))
        else:
            metrics['c1_intensity_weighted_std'] = 0.0

        if metrics['c2_total_intensity'] > 0:
            c2_weighted_dist_sq = (c2_float * distances**2).sum() / metrics['c2_total_intensity']
            metrics['c2_intensity_weighted_std'] = float(np.sqrt(c2_weighted_dist_sq))
        else:
            metrics['c2_intensity_weighted_std'] = 0.0

        # Scatter index: how spread out is the intensity?
        # Higher values = more scattered
        metrics['c1_scatter_index'] = metrics['c1_intensity_weighted_std']
        metrics['c2_scatter_index'] = metrics['c2_intensity_weighted_std']
        metrics['scatter_index_difference'] = abs(
            metrics['c1_scatter_index'] - metrics['c2_scatter_index']
        )

        # 5. RADIAL INTENSITY PROFILES
        # Bin distances and calculate mean intensity per distance bin
        max_dist = distances.max()
        n_bins = 20
        dist_bins = np.linspace(0, max_dist, n_bins + 1)

        c1_radial_profile = []
        c2_radial_profile = []

        for i in range(n_bins):
            bin_mask = (distances >= dist_bins[i]) & (distances < dist_bins[i+1])
            if bin_mask.sum() > 0:
                c1_radial_profile.append(float(c1_float[bin_mask].mean()))
                c2_radial_profile.append(float(c2_float[bin_mask].mean()))
            else:
                c1_radial_profile.append(0.0)
                c2_radial_profile.append(0.0)

        metrics['c1_radial_intensity_profile'] = c1_radial_profile
        metrics['c2_radial_intensity_profile'] = c2_radial_profile

        # Correlate radial profiles
        try:
            profile_r, profile_p = pearsonr(c1_radial_profile, c2_radial_profile)
            metrics['radial_profile_pearson_r'] = float(profile_r)
            metrics['radial_profile_pearson_p'] = float(profile_p)
        except Exception as e:
            metrics['radial_profile_pearson_r'] = None
            metrics['radial_profile_pearson_p'] = None

        # 6. INTENSITY-BASED MANDERS COLOCALIZATION
        # M1 = fraction of C1 intensity colocalized with C2
        # M2 = fraction of C2 intensity colocalized with C1
        # Colocalization = pixels where BOTH channels have intensity above threshold

        # Use Otsu thresholding to define "positive" pixels
        c1_threshold = 0
        c2_threshold = 0

        if c1_float.max() > 0:
            c1_threshold, c1_binary = cv2.threshold(
                c1_img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU
            )

        if c2_float.max() > 0:
            c2_threshold, c2_binary = cv2.threshold(
                c2_img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU
            )

        metrics['c1_otsu_threshold'] = float(c1_threshold)
        metrics['c2_otsu_threshold'] = float(c2_threshold)

        # Define colocalized pixels: both channels above their thresholds
        c1_positive = c1_float > c1_threshold
        c2_positive = c2_float > c2_threshold
        coloc_mask = c1_positive & c2_positive

        # Calculate Manders coefficients
        c1_coloc_intensity = c1_float[coloc_mask].sum() if coloc_mask.sum() > 0 else 0.0
        c2_coloc_intensity = c2_float[coloc_mask].sum() if coloc_mask.sum() > 0 else 0.0

        if metrics['c1_total_intensity'] > 0:
            metrics['manders_m1'] = float(c1_coloc_intensity / metrics['c1_total_intensity'])
        else:
            metrics['manders_m1'] = 0.0

        if metrics['c2_total_intensity'] > 0:
            metrics['manders_m2'] = float(c2_coloc_intensity / metrics['c2_total_intensity'])
        else:
            metrics['manders_m2'] = 0.0

        metrics['colocalized_pixel_count'] = int(coloc_mask.sum())
        metrics['colocalization_fraction'] = float(
            coloc_mask.sum() / (height * width)
        )

        # 7. INTERPRETATION
        # Based on combination of metrics
        pearson_r = metrics.get('pixel_pearson_r', 0)
        manders_avg = (metrics['manders_m1'] + metrics['manders_m2']) / 2
        scatter_diff = metrics['scatter_index_difference']

        # Determine co-scattering interpretation
        if pearson_r is not None and manders_avg is not None:
            if pearson_r > 0.5 and manders_avg > 0.5:
                interpretation = "HIGH_COSCATTERING"
            elif pearson_r > 0.3 and manders_avg > 0.3:
                interpretation = "MODERATE_COSCATTERING"
            elif pearson_r < 0 or manders_avg < 0.2:
                interpretation = "LOW_COSCATTERING"
            else:
                interpretation = "INDETERMINATE"
        else:
            interpretation = "INDETERMINATE"

        metrics['interpretation'] = interpretation
        metrics['status'] = 'SUCCESS'

        return metrics

    def _save_report(self):
        """Save metrics report"""
        report_path = self.output_dir / 'metrics_report.json'

        self.metrics_report['summary'] = {
            'samples_processed': self.stats['samples_processed'],
            'cells_processed': self.stats['cells_processed'],
            'cells_skipped_no_particles': self.stats['cells_skipped']
        }

        with open(report_path, 'w') as f:
            json.dump(self.metrics_report, f, indent=2)

    def _print_results(self):
        """Print metrics results"""
        summary = self.metrics_report['summary']

        click.echo(f"\n✓ Stage 5 completed")
        click.echo(f"  → Output: {self.output_dir}")
        click.echo(f"\n  Samples processed: {summary['samples_processed']}")
        click.echo(f"  Cells processed: {summary['cells_processed']}")
        click.echo(f"  Cells skipped (no particles): {summary['cells_skipped_no_particles']}")

        if len(self.metrics_report['errors']) > 0:
            click.echo(f"\n⚠ {len(self.metrics_report['errors'])} errors. Check metrics_report.json")


def run_stage_5(config, output_dir):
    """Run Stage 5: Scatter Metrics Calculation"""
    calculator = MetricsCalculator(config, output_dir)
    return calculator.run()
