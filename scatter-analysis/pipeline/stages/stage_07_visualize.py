"""
Stage 7: Visualization & Reporting

Purpose:
    Generate visual outputs and comprehensive report.

Input:
    - All previous stage outputs

Output:
    - outputs/07_visualizations/summary_plots/
    - outputs/07_visualizations/report.html
    - outputs/07_visualizations/report.md
"""

import json
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from datetime import datetime
import click


class Visualizer:
    """Creates visualizations and reports"""

    def __init__(self, config, output_dir):
        self.config = config
        self.output_dir = Path(output_dir)
        self.analysis_dir = Path("outputs/06_analysis")

        # Set style
        sns.set_style("whitegrid")
        plt.rcParams['figure.figsize'] = (10, 6)
        plt.rcParams['figure.dpi'] = 100

    def run(self):
        """Execute Stage 7: Visualization & Reporting"""
        click.echo("\n" + "=" * 60)
        click.echo("[Stage 7] Visualization & Reporting")
        click.echo("=" * 60)

        # Create output directories
        plots_dir = self.output_dir / 'summary_plots'
        plots_dir.mkdir(parents=True, exist_ok=True)

        # Load analysis data
        click.echo("\n→ Loading analysis data...")
        stats_file = self.analysis_dir / 'statistics.json'
        with open(stats_file, 'r') as f:
            statistics = json.load(f)

        metrics_file = self.analysis_dir / 'all_metrics.json'
        with open(metrics_file, 'r') as f:
            all_metrics = json.load(f)

        # Generate plots
        click.echo("→ Generating visualizations...")

        # Check if we have temporal data (10 conditions) or binary data (Controls/Treatment)
        groups = list(all_metrics.get('groups', {}).keys())
        if 'Controls' in groups and len(groups) >= 5:
            # Temporal analysis with multiple conditions
            click.echo("   Detected temporal experiment (10 conditions)")
            self._create_temporal_plots(statistics, all_metrics, plots_dir)
        elif 'Controls' in groups and 'Treatment' in groups:
            # Binary comparison (legacy)
            click.echo("   Detected binary comparison (Controls vs Treatment)")
            self._create_comparison_plots(statistics, all_metrics, plots_dir)
        else:
            click.echo(f"   Warning: Unrecognized group structure: {groups}")

        # Generate reports
        click.echo("→ Generating reports...")
        self._create_markdown_report(statistics, all_metrics)
        self._create_html_report(statistics, all_metrics)

        # Print results
        self._print_results()

        return {'plots_dir': str(plots_dir)}

    def _create_temporal_plots(self, statistics, all_metrics, plots_dir):
        """
        Create temporal progression plots for 10-condition experiment

        PhD Thesis Requirement: Show scatter + correlation vs time with error bars
        """

        # Define temporal order (Controls at time=0, then treatments, then recovery)
        condition_order = {
            'Controls': 0,
            '4h treatment': 4,
            '6h treatment': 6,
            '9h treatment': 9,
            '12h treatment': 12,
            '17h treatment': 17,
            '20h treatment': 20,
            'Recovery after 12h': 12.5,  # Slight offset for visualization
            'Recovery after 17h': 17.5,
            'Recovery after 20h': 20.5
        }

        # Extract data by condition
        temporal_data = {}

        for group_name in all_metrics['groups'].keys():
            if group_name not in condition_order:
                continue

            timepoint = condition_order[group_name]
            is_recovery = 'Recovery' in group_name

            # Initialize storage
            temporal_data[timepoint] = {
                'name': group_name,
                'is_recovery': is_recovery,
                'c1_scatter': [],
                'c2_scatter': [],
                'pixel_pearson': [],
                'manders_m1': [],
                'manders_m2': []
            }

            # Extract correlation data from all_metrics
            group_data = all_metrics['groups'][group_name]

            # Load correlation files for this group
            for cell in group_data.get('cells', []):
                if cell.get('status') != 'SUCCESS':
                    continue

                # Get correlation data if available
                if 'correlation' in cell:
                    corr = cell['correlation']
                    temporal_data[timepoint]['pixel_pearson'].append(corr.get('pixel_pearson_r', np.nan))
                    temporal_data[timepoint]['manders_m1'].append(corr.get('manders_m1', np.nan))
                    temporal_data[timepoint]['manders_m2'].append(corr.get('manders_m2', np.nan))
                    temporal_data[timepoint]['c1_scatter'].append(corr.get('c1_scatter_index', np.nan))
                    temporal_data[timepoint]['c2_scatter'].append(corr.get('c2_scatter_index', np.nan))

        # Sort by timepoint
        sorted_timepoints = sorted(temporal_data.keys())

        # Compute means and SEMs
        times = []
        pearson_means = []
        pearson_sems = []
        m1_means = []
        m1_sems = []
        m2_means = []
        m2_sems = []
        c1_scatter_means = []
        c1_scatter_sems = []
        c2_scatter_means = []
        c2_scatter_sems = []
        labels = []
        recovery_mask = []

        for t in sorted_timepoints:
            data = temporal_data[t]
            times.append(t)
            labels.append(data['name'])
            recovery_mask.append(data['is_recovery'])

            # Pearson correlation
            pearson_vals = [v for v in data['pixel_pearson'] if not np.isnan(v)]
            if len(pearson_vals) > 0:
                pearson_means.append(np.mean(pearson_vals))
                pearson_sems.append(np.std(pearson_vals) / np.sqrt(len(pearson_vals)))
            else:
                pearson_means.append(np.nan)
                pearson_sems.append(np.nan)

            # Manders M1
            m1_vals = [v for v in data['manders_m1'] if not np.isnan(v)]
            if len(m1_vals) > 0:
                m1_means.append(np.mean(m1_vals))
                m1_sems.append(np.std(m1_vals) / np.sqrt(len(m1_vals)))
            else:
                m1_means.append(np.nan)
                m1_sems.append(np.nan)

            # Manders M2
            m2_vals = [v for v in data['manders_m2'] if not np.isnan(v)]
            if len(m2_vals) > 0:
                m2_means.append(np.mean(m2_vals))
                m2_sems.append(np.std(m2_vals) / np.sqrt(len(m2_vals)))
            else:
                m2_means.append(np.nan)
                m2_sems.append(np.nan)

            # C1 scatter
            c1_vals = [v for v in data['c1_scatter'] if not np.isnan(v)]
            if len(c1_vals) > 0:
                c1_scatter_means.append(np.mean(c1_vals))
                c1_scatter_sems.append(np.std(c1_vals) / np.sqrt(len(c1_vals)))
            else:
                c1_scatter_means.append(np.nan)
                c1_scatter_sems.append(np.nan)

            # C2 scatter
            c2_vals = [v for v in data['c2_scatter'] if not np.isnan(v)]
            if len(c2_vals) > 0:
                c2_scatter_means.append(np.mean(c2_vals))
                c2_scatter_sems.append(np.std(c2_vals) / np.sqrt(len(c2_vals)))
            else:
                c2_scatter_means.append(np.nan)
                c2_scatter_sems.append(np.nan)

        # Convert to numpy arrays
        times = np.array(times)
        pearson_means = np.array(pearson_means)
        pearson_sems = np.array(pearson_sems)
        m1_means = np.array(m1_means)
        m1_sems = np.array(m1_sems)
        m2_means = np.array(m2_means)
        m2_sems = np.array(m2_sems)
        c1_scatter_means = np.array(c1_scatter_means)
        c1_scatter_sems = np.array(c1_scatter_sems)
        c2_scatter_means = np.array(c2_scatter_means)
        c2_scatter_sems = np.array(c2_scatter_sems)
        recovery_mask = np.array(recovery_mask)

        # --- PLOT 1: Pixel Pearson Correlation vs Time ---
        fig, ax = plt.subplots(figsize=(12, 6))

        # Separate treatment and recovery for color coding
        treatment_idx = ~recovery_mask
        recovery_idx = recovery_mask

        # Plot treatment points
        ax.errorbar(times[treatment_idx], pearson_means[treatment_idx],
                    yerr=pearson_sems[treatment_idx],
                    fmt='o-', color='#e74c3c', linewidth=2, markersize=8,
                    capsize=5, capthick=2, label='Treatment', zorder=3)

        # Plot recovery points
        if np.any(recovery_idx):
            ax.errorbar(times[recovery_idx], pearson_means[recovery_idx],
                       yerr=pearson_sems[recovery_idx],
                       fmt='s-', color='#27ae60', linewidth=2, markersize=8,
                       capsize=5, capthick=2, label='Recovery', zorder=3)

        # Formatting
        ax.set_xlabel('Time (hours)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Pixel Pearson Correlation (C1-C2)', fontsize=12, fontweight='bold')
        ax.set_title('Temporal Progression: Protein X - Structure Y Correlation',
                    fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=11, loc='best')
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

        plt.tight_layout()
        plt.savefig(plots_dir / 'temporal_correlation.png', dpi=300, bbox_inches='tight')
        plt.close()

        # --- PLOT 2: Manders Colocalization vs Time ---
        fig, ax = plt.subplots(figsize=(12, 6))

        # Treatment
        ax.errorbar(times[treatment_idx], m1_means[treatment_idx],
                   yerr=m1_sems[treatment_idx],
                   fmt='o-', color='#3498db', linewidth=2, markersize=8,
                   capsize=5, capthick=2, label='Manders M1 (Treatment)', zorder=3)

        ax.errorbar(times[treatment_idx], m2_means[treatment_idx],
                   yerr=m2_sems[treatment_idx],
                   fmt='^-', color='#9b59b6', linewidth=2, markersize=8,
                   capsize=5, capthick=2, label='Manders M2 (Treatment)', zorder=3)

        # Recovery
        if np.any(recovery_idx):
            ax.errorbar(times[recovery_idx], m1_means[recovery_idx],
                       yerr=m1_sems[recovery_idx],
                       fmt='s-', color='#1abc9c', linewidth=2, markersize=8,
                       capsize=5, capthick=2, label='Manders M1 (Recovery)', zorder=3)

            ax.errorbar(times[recovery_idx], m2_means[recovery_idx],
                       yerr=m2_sems[recovery_idx],
                       fmt='D-', color='#16a085', linewidth=2, markersize=8,
                       capsize=5, capthick=2, label='Manders M2 (Recovery)', zorder=3)

        ax.set_xlabel('Time (hours)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Manders Colocalization Coefficient', fontsize=12, fontweight='bold')
        ax.set_title('Temporal Progression: Colocalization (M1 = C1 overlap with C2, M2 = C2 overlap with C1)',
                    fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10, loc='best')
        ax.set_ylim(0, 1)

        plt.tight_layout()
        plt.savefig(plots_dir / 'temporal_manders.png', dpi=300, bbox_inches='tight')
        plt.close()

        # --- PLOT 3: Scatter Indices vs Time ---
        fig, ax = plt.subplots(figsize=(12, 6))

        # Treatment
        ax.errorbar(times[treatment_idx], c1_scatter_means[treatment_idx],
                   yerr=c1_scatter_sems[treatment_idx],
                   fmt='o-', color='#e67e22', linewidth=2, markersize=8,
                   capsize=5, capthick=2, label='C1 Scatter (Treatment)', zorder=3)

        ax.errorbar(times[treatment_idx], c2_scatter_means[treatment_idx],
                   yerr=c2_scatter_sems[treatment_idx],
                   fmt='^-', color='#d35400', linewidth=2, markersize=8,
                   capsize=5, capthick=2, label='C2 Scatter (Treatment)', zorder=3)

        # Recovery
        if np.any(recovery_idx):
            ax.errorbar(times[recovery_idx], c1_scatter_means[recovery_idx],
                       yerr=c1_scatter_sems[recovery_idx],
                       fmt='s-', color='#2ecc71', linewidth=2, markersize=8,
                       capsize=5, capthick=2, label='C1 Scatter (Recovery)', zorder=3)

            ax.errorbar(times[recovery_idx], c2_scatter_means[recovery_idx],
                       yerr=c2_scatter_sems[recovery_idx],
                       fmt='D-', color='#27ae60', linewidth=2, markersize=8,
                       capsize=5, capthick=2, label='C2 Scatter (Recovery)', zorder=3)

        ax.set_xlabel('Time (hours)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Scatter Index (intensity-weighted std from center)', fontsize=12, fontweight='bold')
        ax.set_title('Temporal Progression: Spatial Scatter', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10, loc='best')

        plt.tight_layout()
        plt.savefig(plots_dir / 'temporal_scatter.png', dpi=300, bbox_inches='tight')
        plt.close()

        click.echo(f"   Created temporal_correlation.png")
        click.echo(f"   Created temporal_manders.png")
        click.echo(f"   Created temporal_scatter.png")

    def _create_comparison_plots(self, statistics, all_metrics, plots_dir):
        """Create comparison plots"""

        # Extract data for plotting
        control_data = {'C1': [], 'C2': [], 'C3': []}
        treatment_data = {'C1': [], 'C2': [], 'C3': []}

        for cell in all_metrics['groups']['Controls']['cells']:
            if cell.get('status') == 'SUCCESS':
                channel = cell.get('channel')
                if channel:
                    control_data[channel].append(cell.get('median_distance_from_center', 0))

        for cell in all_metrics['groups']['Treatment']['cells']:
            if cell.get('status') == 'SUCCESS':
                channel = cell.get('channel')
                if channel:
                    treatment_data[channel].append(cell.get('median_distance_from_center', 0))

        # Plot 1: Box plots by channel
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        fig.suptitle('Median Distance from Center: Controls vs Treatment', fontsize=14, fontweight='bold')

        for idx, channel in enumerate(['C1', 'C2', 'C3']):
            ax = axes[idx]

            data_to_plot = [control_data[channel], treatment_data[channel]]
            bp = ax.boxplot(data_to_plot, labels=['Controls', 'Treatment'],
                            patch_artist=True, widths=0.6)

            # Color boxes
            bp['boxes'][0].set_facecolor('lightblue')
            bp['boxes'][1].set_facecolor('lightcoral')

            ax.set_ylabel('Median Distance (pixels)')
            ax.set_title(f'Channel {channel}')
            ax.grid(True, alpha=0.3)

            # Add sample sizes
            ax.text(1, ax.get_ylim()[1] * 0.95, f'n={len(control_data[channel])}',
                    ha='center', fontsize=9)
            ax.text(2, ax.get_ylim()[1] * 0.95, f'n={len(treatment_data[channel])}',
                    ha='center', fontsize=9)

        plt.tight_layout()
        plt.savefig(plots_dir / 'comparison_boxplots.png', dpi=300, bbox_inches='tight')
        plt.close()

        # Plot 2: Bar chart with error bars
        fig, ax = plt.subplots(figsize=(10, 6))

        channels = ['C1', 'C2', 'C3']
        x = np.arange(len(channels))
        width = 0.35

        control_means = [np.mean(control_data[ch]) if len(control_data[ch]) > 0 else 0 for ch in channels]
        control_stds = [np.std(control_data[ch]) if len(control_data[ch]) > 0 else 0 for ch in channels]
        treatment_means = [np.mean(treatment_data[ch]) if len(treatment_data[ch]) > 0 else 0 for ch in channels]
        treatment_stds = [np.std(treatment_data[ch]) if len(treatment_data[ch]) > 0 else 0 for ch in channels]

        ax.bar(x - width/2, control_means, width, yerr=control_stds, label='Controls',
               color='lightblue', capsize=5)
        ax.bar(x + width/2, treatment_means, width, yerr=treatment_stds, label='Treatment',
               color='lightcoral', capsize=5)

        ax.set_xlabel('Channel')
        ax.set_ylabel('Mean Median Distance (pixels)')
        ax.set_title('Mean Median Distance: Controls vs Treatment', fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(channels)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        plt.savefig(plots_dir / 'comparison_barplot.png', dpi=300, bbox_inches='tight')
        plt.close()

        # Plot 3: Violin plots
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        fig.suptitle('Distribution of Median Distances', fontsize=14, fontweight='bold')

        for idx, channel in enumerate(['C1', 'C2', 'C3']):
            ax = axes[idx]

            data_combined = control_data[channel] + treatment_data[channel]
            groups = ['Controls'] * len(control_data[channel]) + ['Treatment'] * len(treatment_data[channel])

            if len(data_combined) > 0:
                parts = ax.violinplot([control_data[channel], treatment_data[channel]],
                                      positions=[1, 2], widths=0.7, showmeans=True, showmedians=True)

                for pc in parts['bodies']:
                    pc.set_alpha(0.7)

            ax.set_xticks([1, 2])
            ax.set_xticklabels(['Controls', 'Treatment'])
            ax.set_ylabel('Median Distance (pixels)')
            ax.set_title(f'Channel {channel}')
            ax.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        plt.savefig(plots_dir / 'comparison_violinplots.png', dpi=300, bbox_inches='tight')
        plt.close()

    def _create_markdown_report(self, statistics, all_metrics):
        """Create markdown report"""

        report_path = self.output_dir / 'report.md'

        with open(report_path, 'w') as f:
            f.write("# Cell Scatter Analysis Report\n\n")
            f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write("---\n\n")

            # Summary
            f.write("## Executive Summary\n\n")

            # Count cells by group
            groups = all_metrics.get('groups', {})
            total_cells = 0
            for group_name, group_data in groups.items():
                cells_in_group = len([c for c in group_data.get('cells', [])
                                     if c.get('status') == 'SUCCESS'])
                f.write(f"  - **{group_name}**: {cells_in_group} cells\n")
                total_cells += cells_in_group

            f.write(f"\n- **Total Cells Analyzed:** {total_cells}\n\n")

            # Statistical Comparisons
            f.write("## Statistical Comparisons\n\n")
            f.write("See `summary.md` for detailed temporal analysis.\n\n")
            f.write("Key findings from intensity-based correlation analysis:\n\n")

            # Extract correlation stats by group
            if 'by_group' in statistics:
                f.write("### Correlation Metrics by Condition\n\n")
                f.write("| Condition | Pixel Pearson r | Manders M1 | Manders M2 | Sample Size |\n")
                f.write("|-----------|----------------|------------|------------|-------------|\n")

                for group_name in sorted(statistics['by_group'].keys()):
                    group_stats = statistics['by_group'][group_name]
                    if 'correlation' in group_stats:
                        corr = group_stats['correlation']
                        pixel_r = corr.get('pixel_pearson_r', {})
                        m1 = corr.get('manders_m1', {})
                        m2 = corr.get('manders_m2', {})
                        n = pixel_r.get('n', 0)

                        f.write(f"| {group_name} | {pixel_r.get('mean', 0):.3f} ± {pixel_r.get('std', 0):.3f} | {m1.get('mean', 0):.3f} | {m2.get('mean', 0):.3f} | {n} |\n")

                f.write("\n")

            # Interpretation
            f.write("## Interpretation\n\n")
            f.write("Analysis uses intensity-based pixel correlation, not particle positions.\n\n")
            f.write("See visualizations:\n")
            f.write("- `temporal_correlation.png` - Pearson correlation vs time\n")
            f.write("- `temporal_manders.png` - Colocalization coefficients vs time\n")
            f.write("- `temporal_scatter.png` - Scatter indices vs time\n\n")

            f.write("\n---\n\n")
            f.write("*Report generated by Cell Scatter Analysis Pipeline v2.0 (Temporal Analysis)*\n")

    def _create_html_report(self, statistics, all_metrics):
        """Create simple HTML report"""

        report_path = self.output_dir / 'report.html'

        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Cell Scatter Analysis Report - Temporal Co-Scattering</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; background-color: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: auto; background-color: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        h3 {{ color: #7f8c8d; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        tr:nth-child(even) {{ background-color: #f9f9f9; }}
        .stat {{ display: inline-block; margin: 10px 20px; }}
        .stat-label {{ font-weight: bold; color: #7f8c8d; }}
        .stat-value {{ font-size: 24px; color: #2c3e50; }}
        .significant {{ color: #27ae60; font-weight: bold; }}
        .not-significant {{ color: #e74c3c; }}
        img {{ max-width: 100%; height: auto; margin: 20px 0; border: 1px solid #ddd; border-radius: 4px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Temporal Co-Scattering Analysis Report</h1>
        <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p><strong>Analysis Type:</strong> Intensity-based pixel correlation (10 temporal conditions)</p>

        <h2>Summary</h2>
"""

        # Count cells by group
        groups = all_metrics.get('groups', {})
        total_cells = 0
        for group_name, group_data in groups.items():
            cells_in_group = len([c for c in group_data.get('cells', [])
                                 if c.get('status') == 'SUCCESS'])
            html += f"""
        <div class="stat">
            <div class="stat-label">{group_name}</div>
            <div class="stat-value">{cells_in_group}</div>
        </div>"""
            total_cells += cells_in_group

        html += f"""
        <div class="stat">
            <div class="stat-label"><strong>Total Cells</strong></div>
            <div class="stat-value"><strong>{total_cells}</strong></div>
        </div>

        <h2>Temporal Progression Visualizations</h2>
        <p>Analysis shows how protein X (C2) and structure Y (C1) co-scatter across drug treatment time course.</p>

        <h3>Correlation vs Time</h3>
        <img src="summary_plots/temporal_correlation.png" alt="Temporal Correlation">
        <p>Pixel-based Pearson correlation between C1 and C2 intensities. Higher values indicate more coordinated distribution.</p>

        <h3>Colocalization vs Time</h3>
        <img src="summary_plots/temporal_manders.png" alt="Temporal Manders">
        <p>Manders coefficients: M1 = fraction of C1 overlapping with C2, M2 = fraction of C2 overlapping with C1.</p>

        <h3>Spatial Scatter vs Time</h3>
        <img src="summary_plots/temporal_scatter.png" alt="Temporal Scatter">
        <p>Scatter index measures intensity-weighted standard deviation from cell center. Higher values = more dispersed.</p>

        <h2>Statistical Details</h2>
        <p>See <code>summary.md</code> for complete statistical analysis with Bonferroni correction.</p>
        <p>See <code>outputs/06_analysis/statistics.json</code> for raw statistical data.</p>

        <hr>
        <p><em>Report generated by Cell Scatter Analysis Pipeline v2.0 (Temporal Co-Scattering)</em></p>
    </div>
</body>
</html>
"""

        with open(report_path, 'w') as f:
            f.write(html)

    def _print_results(self):
        """Print results"""

        click.echo(f"\n✓ Stage 7 completed")
        click.echo(f"  → Output: {self.output_dir}")
        click.echo(f"\n  Generated:")
        click.echo(f"    - summary_plots/ (3 plots)")
        click.echo(f"    - report.html")
        click.echo(f"    - report.md")
        click.echo(f"\n  Open report:")
        click.echo(f"    open {self.output_dir / 'report.html'}")


def run_stage_7(config, output_dir):
    """Run Stage 7: Visualization & Reporting"""
    visualizer = Visualizer(config, output_dir)
    return visualizer.run()
