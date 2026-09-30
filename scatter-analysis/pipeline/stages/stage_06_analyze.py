"""
Stage 6: Aggregation & Comparison

Purpose:
    Aggregate all metrics and perform statistical comparisons.

Input:
    - outputs/05_metrics/**/*.json

Output:
    - outputs/06_analysis/all_cells_summary.csv
    - outputs/06_analysis/all_metrics.json
    - outputs/06_analysis/statistics.json
    - outputs/06_analysis/comparison.json
"""

import json
import csv
import numpy as np
from pathlib import Path
from datetime import datetime
from scipy import stats
from tqdm import tqdm
import click


class DataAnalyzer:
    """Aggregates and analyzes scatter metrics"""

    def __init__(self, config, output_dir):
        self.config = config
        self.output_dir = Path(output_dir)
        self.metrics_dir = Path("outputs/05_metrics")

        # Dynamically discover groups (conditions) from metrics directory
        self.groups = self._discover_groups()

        self.all_metrics = {
            'created_at': datetime.now().isoformat(),
            'groups': {group: {'cells': [], 'correlation_cells': []} for group in self.groups}
        }

        self.statistics = {
            'created_at': datetime.now().isoformat(),
            'by_group': {},
            'comparison': {}
        }

    def _discover_groups(self):
        """Dynamically discover condition groups from metrics directory"""
        groups = []
        if self.metrics_dir.exists():
            for item in self.metrics_dir.iterdir():
                if item.is_dir():
                    groups.append(item.name)

        # Fallback if no groups found
        if not groups:
            groups = ['Controls', 'Treatment']

        return sorted(groups)

    def run(self):
        """Execute Stage 6: Aggregation & Comparison"""
        click.echo("\n" + "=" * 60)
        click.echo("[Stage 6] Aggregation & Comparison")
        click.echo("=" * 60)

        # Collect all metrics
        click.echo("\n→ Collecting all cell metrics...")
        self._collect_all_metrics()

        # Calculate group statistics
        click.echo("→ Calculating group statistics...")
        self._calculate_group_statistics()

        # Perform comparisons
        click.echo("→ Performing statistical comparisons...")
        self._perform_comparisons()

        # Apply multiple comparison correction
        self._apply_multiple_comparison_correction()

        # Save outputs
        click.echo("→ Saving aggregated data...")
        self._save_outputs()

        # Print results
        self._print_results()

        return self.all_metrics, self.statistics

    def _collect_all_metrics(self):
        """Collect all metrics from Stage 5 (per-channel + correlation)"""

        for group in self.groups:
            group_dir = self.metrics_dir / group

            if not group_dir.exists():
                continue

            # Iterate through samples
            for sample_dir in group_dir.iterdir():
                if not sample_dir.is_dir():
                    continue

                sample_name = sample_dir.name

                # Collect per-channel metrics
                for channel in ['C1', 'C2', 'C3']:
                    channel_dir = sample_dir / channel

                    if not channel_dir.exists():
                        continue

                    # Iterate through cells
                    for metrics_file in channel_dir.glob('cell_*_metrics.json'):
                        with open(metrics_file, 'r') as f:
                            metrics = json.load(f)

                        # Add group and sample info
                        metrics['group'] = group
                        metrics['sample_name'] = sample_name

                        self.all_metrics['groups'][group]['cells'].append(metrics)

                # Collect correlation metrics (NEW!)
                correlation_dir = sample_dir / 'correlation'
                if correlation_dir.exists():
                    for corr_file in correlation_dir.glob('cell_*_correlation.json'):
                        with open(corr_file, 'r') as f:
                            corr_metrics = json.load(f)

                        # Add group and sample info
                        corr_metrics['group'] = group
                        corr_metrics['sample_name'] = sample_name

                        self.all_metrics['groups'][group]['correlation_cells'].append(corr_metrics)

    def _calculate_group_statistics(self):
        """Calculate statistics for each group (including correlation metrics)"""

        for group in self.groups:
            cells = self.all_metrics['groups'][group]['cells']

            # Group by channel
            by_channel = {'C1': [], 'C2': [], 'C3': []}

            for cell in cells:
                if cell.get('status') != 'SUCCESS':
                    continue

                channel = cell.get('channel')
                if channel in by_channel:
                    by_channel[channel].append(cell)

            # Calculate statistics for each channel
            channel_stats = {}

            for channel, channel_cells in by_channel.items():
                if len(channel_cells) == 0:
                    continue

                # Extract key metrics
                median_distances = [c.get('median_distance_from_center', 0) for c in channel_cells]
                mean_distances = [c.get('mean_distance_from_center', 0) for c in channel_cells]
                particle_counts = [c.get('total_particles', 0) for c in channel_cells]
                convex_hull_areas = [c.get('convex_hull_area', 0) for c in channel_cells]

                channel_stats[channel] = {
                    'cell_count': len(channel_cells),
                    'median_distance': {
                        'mean': float(np.mean(median_distances)),
                        'median': float(np.median(median_distances)),
                        'std': float(np.std(median_distances)),
                        'min': float(np.min(median_distances)),
                        'max': float(np.max(median_distances))
                    },
                    'mean_distance': {
                        'mean': float(np.mean(mean_distances)),
                        'median': float(np.median(mean_distances)),
                        'std': float(np.std(mean_distances))
                    },
                    'particle_count': {
                        'mean': float(np.mean(particle_counts)),
                        'median': float(np.median(particle_counts)),
                        'std': float(np.std(particle_counts))
                    },
                    'convex_hull_area': {
                        'mean': float(np.mean(convex_hull_areas)),
                        'median': float(np.median(convex_hull_areas)),
                        'std': float(np.std(convex_hull_areas))
                    }
                }

            # Calculate correlation statistics (NEW!)
            corr_cells = self.all_metrics['groups'][group]['correlation_cells']
            corr_cells_success = [c for c in corr_cells if c.get('status') == 'SUCCESS']

            if len(corr_cells_success) > 0:
                # Use intensity-based metrics (pixel_pearson_r, pixel_spearman_rho)
                pearson_rs = [c.get('pixel_pearson_r') for c in corr_cells_success
                              if c.get('pixel_pearson_r') is not None]
                spearman_rhos = [c.get('pixel_spearman_rho') for c in corr_cells_success
                                  if c.get('pixel_spearman_rho') is not None]
                manders_m1s = [c.get('manders_m1', 0) for c in corr_cells_success]
                manders_m2s = [c.get('manders_m2', 0) for c in corr_cells_success]
                # Use intensity_centroid_distance instead of centroid_distance
                centroid_dists = [c.get('intensity_centroid_distance', 0) for c in corr_cells_success]
                # Also collect scatter indices
                c1_scatter = [c.get('c1_scatter_index', 0) for c in corr_cells_success
                             if c.get('c1_scatter_index') is not None]
                c2_scatter = [c.get('c2_scatter_index', 0) for c in corr_cells_success
                             if c.get('c2_scatter_index') is not None]

                channel_stats['correlation'] = {
                    'cell_count': len(corr_cells_success),
                    'pearson_r': {
                        'mean': float(np.mean(pearson_rs)) if pearson_rs else None,
                        'median': float(np.median(pearson_rs)) if pearson_rs else None,
                        'std': float(np.std(pearson_rs)) if pearson_rs else None
                    },
                    'spearman_rho': {
                        'mean': float(np.mean(spearman_rhos)) if spearman_rhos else None,
                        'median': float(np.median(spearman_rhos)) if spearman_rhos else None,
                        'std': float(np.std(spearman_rhos)) if spearman_rhos else None
                    },
                    'manders_m1': {
                        'mean': float(np.mean(manders_m1s)),
                        'median': float(np.median(manders_m1s)),
                        'std': float(np.std(manders_m1s))
                    },
                    'manders_m2': {
                        'mean': float(np.mean(manders_m2s)),
                        'median': float(np.median(manders_m2s)),
                        'std': float(np.std(manders_m2s))
                    },
                    'centroid_distance': {
                        'mean': float(np.mean(centroid_dists)),
                        'median': float(np.median(centroid_dists)),
                        'std': float(np.std(centroid_dists))
                    },
                    'c1_scatter_index': {
                        'mean': float(np.mean(c1_scatter)) if c1_scatter else None,
                        'median': float(np.median(c1_scatter)) if c1_scatter else None,
                        'std': float(np.std(c1_scatter)) if c1_scatter else None
                    },
                    'c2_scatter_index': {
                        'mean': float(np.mean(c2_scatter)) if c2_scatter else None,
                        'median': float(np.median(c2_scatter)) if c2_scatter else None,
                        'std': float(np.std(c2_scatter)) if c2_scatter else None
                    }
                }

            self.statistics['by_group'][group] = channel_stats

    def _perform_comparisons(self):
        """Perform statistical comparisons (all conditions vs Controls)"""

        # Determine control group (should contain 'Control' or 'Untreated' in name)
        control_group = None
        for group in self.groups:
            if 'control' in group.lower() or 'untreated' in group.lower():
                control_group = group
                break

        if not control_group:
            click.echo("⚠ Warning: No control group found. Skipping comparisons.")
            return

        # Initialize comparison structure
        # Format: comparison[group][channel] or comparison[group]['correlation']
        self.statistics['comparison'] = {}

        # Compare each non-control group to control
        for group in self.groups:
            if group == control_group:
                continue  # Skip control vs control

            self.statistics['comparison'][group] = {}

            # Compare per-channel metrics (C1, C2, C3)
            for channel in ['C1', 'C2', 'C3']:
                # Get control cells for this channel
                control_cells = [c for c in self.all_metrics['groups'][control_group]['cells']
                                 if c.get('channel') == channel and c.get('status') == 'SUCCESS']

                # Get treatment cells for this channel
                treatment_cells = [c for c in self.all_metrics['groups'][group]['cells']
                                   if c.get('channel') == channel and c.get('status') == 'SUCCESS']

                if len(control_cells) == 0 or len(treatment_cells) == 0:
                    continue

                # Extract median distances
                control_medians = [c.get('median_distance_from_center', 0) for c in control_cells]
                treatment_medians = [c.get('median_distance_from_center', 0) for c in treatment_cells]

                # Perform t-test
                t_stat, p_value = stats.ttest_ind(control_medians, treatment_medians)

                # Calculate effect size (Cohen's d)
                pooled_std = np.sqrt((np.std(control_medians)**2 + np.std(treatment_medians)**2) / 2)
                if pooled_std > 0:
                    cohens_d = (np.mean(treatment_medians) - np.mean(control_medians)) / pooled_std
                else:
                    cohens_d = 0.0

                # Calculate difference
                control_mean = np.mean(control_medians)
                treatment_mean = np.mean(treatment_medians)
                difference = treatment_mean - control_mean
                percent_change = (difference / control_mean * 100) if control_mean > 0 else 0.0

                self.statistics['comparison'][group][channel] = {
                    'control_mean': float(control_mean),
                    'treatment_mean': float(treatment_mean),
                    'difference': float(difference),
                    'percent_change': float(percent_change),
                    't_statistic': float(t_stat),
                    'p_value': float(p_value),
                    'cohens_d': float(cohens_d),
                    'significance': 'significant' if p_value < 0.05 else 'not_significant',
                    'control_n': len(control_cells),
                    'treatment_n': len(treatment_cells)
                }

            # Compare correlation metrics (C1-C2 co-scattering: Protein X vs Structure Y)
            control_corr_cells = [c for c in self.all_metrics['groups'][control_group].get('correlation_cells', [])
                                  if c.get('status') == 'SUCCESS']
            treatment_corr_cells = [c for c in self.all_metrics['groups'][group].get('correlation_cells', [])
                                    if c.get('status') == 'SUCCESS']

            if len(control_corr_cells) > 0 and len(treatment_corr_cells) > 0:
                # Extract intensity-based correlation metrics
                control_pearson = [c.get('pixel_pearson_r') for c in control_corr_cells
                                   if c.get('pixel_pearson_r') is not None]
                treatment_pearson = [c.get('pixel_pearson_r') for c in treatment_corr_cells
                                     if c.get('pixel_pearson_r') is not None]

                control_manders_m1 = [c.get('manders_m1', 0) for c in control_corr_cells]
                treatment_manders_m1 = [c.get('manders_m1', 0) for c in treatment_corr_cells]

                control_manders_m2 = [c.get('manders_m2', 0) for c in control_corr_cells]
                treatment_manders_m2 = [c.get('manders_m2', 0) for c in treatment_corr_cells]

                # Extract scatter indices
                control_c1_scatter = [c.get('c1_scatter_index') for c in control_corr_cells
                                     if c.get('c1_scatter_index') is not None]
                treatment_c1_scatter = [c.get('c1_scatter_index') for c in treatment_corr_cells
                                       if c.get('c1_scatter_index') is not None]

                control_c2_scatter = [c.get('c2_scatter_index') for c in control_corr_cells
                                     if c.get('c2_scatter_index') is not None]
                treatment_c2_scatter = [c.get('c2_scatter_index') for c in treatment_corr_cells
                                       if c.get('c2_scatter_index') is not None]

                control_centroid_dist = [c.get('intensity_centroid_distance', 0) for c in control_corr_cells]
                treatment_centroid_dist = [c.get('intensity_centroid_distance', 0) for c in treatment_corr_cells]

                # Statistical tests for correlation metrics
                correlation_comparison = {}

                # Pearson correlation comparison
                if len(control_pearson) > 0 and len(treatment_pearson) > 0:
                    t_stat, p_value = stats.ttest_ind(control_pearson, treatment_pearson)
                    correlation_comparison['pearson_r'] = {
                        'control_mean': float(np.mean(control_pearson)),
                        'treatment_mean': float(np.mean(treatment_pearson)),
                        'control_std': float(np.std(control_pearson)),
                        'treatment_std': float(np.std(treatment_pearson)),
                        'difference': float(np.mean(treatment_pearson) - np.mean(control_pearson)),
                        't_statistic': float(t_stat),
                        'p_value': float(p_value),
                        'significance': 'significant' if p_value < 0.05 else 'not_significant',
                        'control_n': len(control_pearson),
                        'treatment_n': len(treatment_pearson)
                    }

                # Manders M1 comparison
                if len(control_manders_m1) > 0 and len(treatment_manders_m1) > 0:
                    t_stat, p_value = stats.ttest_ind(control_manders_m1, treatment_manders_m1)
                    correlation_comparison['manders_m1'] = {
                        'control_mean': float(np.mean(control_manders_m1)),
                        'treatment_mean': float(np.mean(treatment_manders_m1)),
                        'control_std': float(np.std(control_manders_m1)),
                        'treatment_std': float(np.std(treatment_manders_m1)),
                        'difference': float(np.mean(treatment_manders_m1) - np.mean(control_manders_m1)),
                        't_statistic': float(t_stat),
                        'p_value': float(p_value),
                        'significance': 'significant' if p_value < 0.05 else 'not_significant',
                        'control_n': len(control_manders_m1),
                        'treatment_n': len(treatment_manders_m1)
                    }

                # Manders M2 comparison
                if len(control_manders_m2) > 0 and len(treatment_manders_m2) > 0:
                    t_stat, p_value = stats.ttest_ind(control_manders_m2, treatment_manders_m2)
                    correlation_comparison['manders_m2'] = {
                        'control_mean': float(np.mean(control_manders_m2)),
                        'treatment_mean': float(np.mean(treatment_manders_m2)),
                        'control_std': float(np.std(control_manders_m2)),
                        'treatment_std': float(np.std(treatment_manders_m2)),
                        'difference': float(np.mean(treatment_manders_m2) - np.mean(control_manders_m2)),
                        't_statistic': float(t_stat),
                        'p_value': float(p_value),
                        'significance': 'significant' if p_value < 0.05 else 'not_significant',
                        'control_n': len(control_manders_m2),
                        'treatment_n': len(treatment_manders_m2)
                    }

                # Centroid distance comparison
                if len(control_centroid_dist) > 0 and len(treatment_centroid_dist) > 0:
                    t_stat, p_value = stats.ttest_ind(control_centroid_dist, treatment_centroid_dist)
                    correlation_comparison['centroid_distance'] = {
                        'control_mean': float(np.mean(control_centroid_dist)),
                        'treatment_mean': float(np.mean(treatment_centroid_dist)),
                        'control_std': float(np.std(control_centroid_dist)),
                        'treatment_std': float(np.std(treatment_centroid_dist)),
                        'difference': float(np.mean(treatment_centroid_dist) - np.mean(control_centroid_dist)),
                        't_statistic': float(t_stat),
                        'p_value': float(p_value),
                        'significance': 'significant' if p_value < 0.05 else 'not_significant',
                        'control_n': len(control_centroid_dist),
                        'treatment_n': len(treatment_centroid_dist)
                    }

                # C1 scatter index comparison (SCATTER HYPOTHESIS TEST)
                if len(control_c1_scatter) > 0 and len(treatment_c1_scatter) > 0:
                    t_stat, p_value = stats.ttest_ind(control_c1_scatter, treatment_c1_scatter)

                    # Calculate Cohen's d effect size
                    pooled_std = np.sqrt(((len(control_c1_scatter) - 1) * np.std(control_c1_scatter)**2 +
                                          (len(treatment_c1_scatter) - 1) * np.std(treatment_c1_scatter)**2) /
                                         (len(control_c1_scatter) + len(treatment_c1_scatter) - 2))
                    cohens_d = (np.mean(treatment_c1_scatter) - np.mean(control_c1_scatter)) / pooled_std if pooled_std > 0 else 0

                    correlation_comparison['c1_scatter_index'] = {
                        'control_mean': float(np.mean(control_c1_scatter)),
                        'treatment_mean': float(np.mean(treatment_c1_scatter)),
                        'control_std': float(np.std(control_c1_scatter)),
                        'treatment_std': float(np.std(treatment_c1_scatter)),
                        'difference': float(np.mean(treatment_c1_scatter) - np.mean(control_c1_scatter)),
                        't_statistic': float(t_stat),
                        'p_value': float(p_value),
                        'cohens_d': float(cohens_d),
                        'significance': 'significant' if p_value < 0.05 else 'not_significant',
                        'control_n': len(control_c1_scatter),
                        'treatment_n': len(treatment_c1_scatter)
                    }

                # C2 scatter index comparison (SCATTER HYPOTHESIS TEST - PRIMARY METRIC)
                if len(control_c2_scatter) > 0 and len(treatment_c2_scatter) > 0:
                    t_stat, p_value = stats.ttest_ind(control_c2_scatter, treatment_c2_scatter)

                    # Calculate Cohen's d effect size
                    pooled_std = np.sqrt(((len(control_c2_scatter) - 1) * np.std(control_c2_scatter)**2 +
                                          (len(treatment_c2_scatter) - 1) * np.std(treatment_c2_scatter)**2) /
                                         (len(control_c2_scatter) + len(treatment_c2_scatter) - 2))
                    cohens_d = (np.mean(treatment_c2_scatter) - np.mean(control_c2_scatter)) / pooled_std if pooled_std > 0 else 0

                    correlation_comparison['c2_scatter_index'] = {
                        'control_mean': float(np.mean(control_c2_scatter)),
                        'treatment_mean': float(np.mean(treatment_c2_scatter)),
                        'control_std': float(np.std(control_c2_scatter)),
                        'treatment_std': float(np.std(treatment_c2_scatter)),
                        'difference': float(np.mean(treatment_c2_scatter) - np.mean(control_c2_scatter)),
                        't_statistic': float(t_stat),
                        'p_value': float(p_value),
                        'cohens_d': float(cohens_d),
                        'significance': 'significant' if p_value < 0.05 else 'not_significant',
                        'control_n': len(control_c2_scatter),
                        'treatment_n': len(treatment_c2_scatter)
                    }

                self.statistics['comparison'][group]['correlation'] = correlation_comparison

    def _apply_multiple_comparison_correction(self):
        """Apply Bonferroni correction for multiple comparisons

        With 9 experimental conditions compared to Controls, we need to correct
        for multiple comparisons to maintain PhD-level statistical rigor.

        Bonferroni correction: p_corrected = p_raw * n_comparisons
        Conservative but appropriate for thesis work.
        """

        if not self.statistics.get('comparison'):
            return

        # Count total number of statistical tests performed
        n_comparisons = 0
        all_p_values = []

        # Collect all p-values from channel comparisons
        for group, group_comparisons in self.statistics['comparison'].items():
            for channel in ['C1', 'C2', 'C3']:
                if channel in group_comparisons:
                    n_comparisons += 1
                    all_p_values.append(group_comparisons[channel]['p_value'])

            # Also count correlation comparisons
            if 'correlation' in group_comparisons:
                corr_comp = group_comparisons['correlation']
                for metric in ['pearson_r', 'manders_m1', 'manders_m2', 'centroid_distance',
                              'c1_scatter_index', 'c2_scatter_index']:
                    if metric in corr_comp:
                        n_comparisons += 1
                        all_p_values.append(corr_comp[metric]['p_value'])

        if n_comparisons == 0:
            return

        # Apply Bonferroni correction
        click.echo(f"\n→ Applying Bonferroni correction for {n_comparisons} comparisons...")

        for group, group_comparisons in self.statistics['comparison'].items():
            # Correct channel comparisons
            for channel in ['C1', 'C2', 'C3']:
                if channel in group_comparisons:
                    comp = group_comparisons[channel]
                    p_raw = comp['p_value']
                    p_corrected = min(p_raw * n_comparisons, 1.0)  # Cap at 1.0

                    comp['p_value_raw'] = p_raw
                    comp['p_value_corrected'] = p_corrected
                    comp['p_value'] = p_corrected  # Update main p_value
                    comp['significance_raw'] = 'significant' if p_raw < 0.05 else 'not_significant'
                    comp['significance'] = 'significant' if p_corrected < 0.05 else 'not_significant'
                    comp['multiple_comparison_correction'] = 'bonferroni'
                    comp['n_comparisons'] = n_comparisons

            # Correct correlation comparisons
            if 'correlation' in group_comparisons:
                corr_comp = group_comparisons['correlation']
                for metric in ['pearson_r', 'manders_m1', 'manders_m2', 'centroid_distance',
                              'c1_scatter_index', 'c2_scatter_index']:
                    if metric in corr_comp:
                        p_raw = corr_comp[metric]['p_value']
                        p_corrected = min(p_raw * n_comparisons, 1.0)

                        corr_comp[metric]['p_value_raw'] = p_raw
                        corr_comp[metric]['p_value_corrected'] = p_corrected
                        corr_comp[metric]['p_value'] = p_corrected
                        corr_comp[metric]['significance_raw'] = 'significant' if p_raw < 0.05 else 'not_significant'
                        corr_comp[metric]['significance'] = 'significant' if p_corrected < 0.05 else 'not_significant'
                        corr_comp[metric]['multiple_comparison_correction'] = 'bonferroni'
                        corr_comp[metric]['n_comparisons'] = n_comparisons

        # Store correction metadata
        self.statistics['multiple_comparison_correction'] = {
            'method': 'bonferroni',
            'n_comparisons': n_comparisons,
            'alpha': 0.05,
            'corrected_alpha': 0.05 / n_comparisons,
            'rationale': 'Conservative correction appropriate for PhD thesis with multiple temporal conditions'
        }

    def _save_outputs(self):
        """Save all output files"""

        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 1. Save all_metrics.json
        metrics_file = self.output_dir / 'all_metrics.json'
        with open(metrics_file, 'w') as f:
            json.dump(self.all_metrics, f, indent=2)

        # 2. Save statistics.json
        stats_file = self.output_dir / 'statistics.json'
        with open(stats_file, 'w') as f:
            json.dump(self.statistics, f, indent=2)

        # 3. Save all_cells_summary.csv (flat table)
        csv_file = self.output_dir / 'all_cells_summary.csv'
        self._save_csv(csv_file)

        # 4. Save comparison.json
        comparison_file = self.output_dir / 'comparison.json'
        with open(comparison_file, 'w') as f:
            json.dump(self.statistics['comparison'], f, indent=2)

    def _save_csv(self, csv_file):
        """Save flat CSV file with all cells"""

        fieldnames = [
            'group', 'sample_name', 'cell_id', 'channel',
            'total_particles', 'median_distance_from_center', 'mean_distance_from_center',
            'std_distance_from_center', 'max_distance_from_center',
            'convex_hull_area', 'particle_density', 'percent_cell_occupied',
            'outliers_detected', 'metric_used', 'status'
        ]

        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
            writer.writeheader()

            # Iterate through all discovered groups
            for group in self.groups:
                for cell in self.all_metrics['groups'][group]['cells']:
                    writer.writerow(cell)

    def _print_results(self):
        """Print analysis results"""

        click.echo(f"\n✓ Stage 6 completed")
        click.echo(f"  → Output: {self.output_dir}")

        # Print group summaries
        click.echo(f"\n  Group Statistics:")
        for group, group_stats in self.statistics['by_group'].items():
            # Calculate total cells (excluding correlation which is not a channel)
            total_cells = sum(ch.get('cell_count', 0) for ch_name, ch in group_stats.items()
                              if ch_name in ['C1', 'C2', 'C3'])
            click.echo(f"\n  {group}: {total_cells} cells analyzed")

            for channel_name, ch_stats in group_stats.items():
                if channel_name == 'correlation':
                    # Print correlation stats
                    click.echo(f"    Correlation (C1-C2 Co-scattering): {ch_stats['cell_count']} cells")
                    if ch_stats.get('pearson_r', {}).get('mean') is not None:
                        pearson_mean = ch_stats['pearson_r']['mean']
                        pearson_std = ch_stats['pearson_r']['std']
                        click.echo(f"      Pearson r: {pearson_mean:.3f} ± {pearson_std:.3f}")
                    if ch_stats.get('manders_m1', {}).get('mean') is not None:
                        m1_mean = ch_stats['manders_m1']['mean']
                        m2_mean = ch_stats['manders_m2']['mean']
                        click.echo(f"      Manders M1/M2: {m1_mean:.3f} / {m2_mean:.3f}")
                else:
                    # Print per-channel stats
                    median_dist = ch_stats['median_distance']['mean']
                    click.echo(f"    {channel_name}: {ch_stats['cell_count']} cells, "
                               f"median distance = {median_dist:.2f} ± {ch_stats['median_distance']['std']:.2f}")

        # Print comparisons
        if len(self.statistics.get('comparison', {})) > 0:
            click.echo(f"\n  Statistical Comparisons (All Conditions vs Controls):")

            for group, group_comparisons in self.statistics['comparison'].items():
                click.echo(f"\n  {group}:")

                # Print per-channel comparisons
                for channel in ['C1', 'C2', 'C3']:
                    if channel in group_comparisons:
                        comp = group_comparisons[channel]
                        click.echo(f"\n    {channel}:")
                        click.echo(f"      Controls:  {comp['control_mean']:.2f} (n={comp['control_n']})")
                        click.echo(f"      {group}: {comp['treatment_mean']:.2f} (n={comp['treatment_n']})")
                        click.echo(f"      Difference: {comp['difference']:.2f} ({comp['percent_change']:.1f}%)")
                        click.echo(f"      p-value: {comp['p_value']:.6f} ({comp['significance']})")
                        click.echo(f"      Effect size (Cohen's d): {comp['cohens_d']:.2f}")

                # Print correlation comparisons
                if 'correlation' in group_comparisons:
                    click.echo(f"\n    Correlation (C1-C2 Co-scattering):")
                    corr_comp = group_comparisons['correlation']

                    if 'pearson_r' in corr_comp:
                        pr = corr_comp['pearson_r']
                        click.echo(f"      Pearson r:")
                        click.echo(f"        Controls: {pr['control_mean']:.3f} ± {pr['control_std']:.3f} (n={pr['control_n']})")
                        click.echo(f"        {group}: {pr['treatment_mean']:.3f} ± {pr['treatment_std']:.3f} (n={pr['treatment_n']})")
                        click.echo(f"        Difference: {pr['difference']:.3f} (p={pr['p_value']:.6f}, {pr['significance']})")

                    if 'manders_m1' in corr_comp:
                        m1 = corr_comp['manders_m1']
                        click.echo(f"      Manders M1 (C2 overlap with C3):")
                        click.echo(f"        Controls: {m1['control_mean']:.3f} ± {m1['control_std']:.3f}")
                        click.echo(f"        {group}: {m1['treatment_mean']:.3f} ± {m1['treatment_std']:.3f}")
                        click.echo(f"        Difference: {m1['difference']:.3f} (p={m1['p_value']:.6f}, {m1['significance']})")

                    if 'centroid_distance' in corr_comp:
                        cd = corr_comp['centroid_distance']
                        click.echo(f"      Centroid distance:")
                        click.echo(f"        Controls: {cd['control_mean']:.2f} ± {cd['control_std']:.2f} pixels")
                        click.echo(f"        {group}: {cd['treatment_mean']:.2f} ± {cd['treatment_std']:.2f} pixels")
                        click.echo(f"        Difference: {cd['difference']:.2f} (p={cd['p_value']:.6f}, {cd['significance']})")

            # Print multiple comparison correction info
            if 'multiple_comparison_correction' in self.statistics:
                mcc = self.statistics['multiple_comparison_correction']
                click.echo(f"\n  Multiple Comparison Correction:")
                click.echo(f"    Method: {mcc['method'].capitalize()}")
                click.echo(f"    Number of comparisons: {mcc['n_comparisons']}")
                click.echo(f"    Corrected alpha: {mcc['corrected_alpha']:.6f} (original: {mcc['alpha']})")
                click.echo(f"    Note: All p-values shown above are Bonferroni-corrected")


def run_stage_6(config, output_dir):
    """Run Stage 6: Aggregation & Comparison"""
    analyzer = DataAnalyzer(config, output_dir)
    return analyzer.run()
