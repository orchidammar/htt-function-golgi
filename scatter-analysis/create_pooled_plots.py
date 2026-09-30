#!/usr/bin/env python3
"""
Create publication-quality plots for pooled experiment analysis
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib.patches import Rectangle
from scipy import stats as scipy_stats

# Set publication style
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("Set2")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 9

def load_pooled_data():
    """Load scatter data from all three experiments"""

    def load_from_experiment(output_dir):
        """Load scatter indices from one experiment"""
        metrics_dir = Path(output_dir) / "05_metrics"
        data = {}

        for condition_dir in metrics_dir.iterdir():
            if not condition_dir.is_dir():
                continue

            condition_name = condition_dir.name
            data[condition_name] = {'c1': [], 'c2': []}

            for sample_dir in condition_dir.iterdir():
                if not sample_dir.is_dir():
                    continue

                correlation_dir = sample_dir / "correlation"
                if not correlation_dir.exists():
                    continue

                for corr_file in correlation_dir.glob("*_correlation.json"):
                    with open(corr_file, 'r') as f:
                        corr_data = json.load(f)

                    if corr_data.get('c1_scatter_index'):
                        data[condition_name]['c1'].append(corr_data['c1_scatter_index'])
                    if corr_data.get('c2_scatter_index'):
                        data[condition_name]['c2'].append(corr_data['c2_scatter_index'])

        return data

    print("Loading data from all experiments...")
    exp1_data = load_from_experiment("experiments/outputs")
    exp2_data = load_from_experiment("experiments/outputs_exp2")
    exp3_data = load_from_experiment("experiments/outputs_exp3")

    return exp1_data, exp2_data, exp3_data

def map_condition_to_timepoint(condition_name):
    """Map condition names to standardized time points"""
    mapping = {
        'Controls': ('Controls', 0),
        '4 hours treatment': ('4h', 4),
        '6 hours treatment': ('6h', 6),
        '9 hours treatment': ('9h', 9),
        '12 hours treatment': ('12h', 12),
        '17 hours treatment': ('17h', 17),
        '20 hours treatment': ('20h', 20),
        'Recovery after 12h treatment': ('Recovery 12h', 12.5),
        'Recovery after 12 hours': ('Recovery 12h', 12.5),
        'Recovery after 17h treatment': ('Recovery 17h', 17.5),
        'Recovery after 17 hours': ('Recovery 17h', 17.5),
        'Recovery after 20h treatment': ('Recovery 20h', 20.5),
        'Recovery after 20 hours': ('Recovery 20h', 20.5),
    }
    return mapping.get(condition_name, (condition_name, None))

def pool_by_timepoint(exp1_data, exp2_data, exp3_data):
    """Pool data across experiments by time point"""
    pooled = {}

    for exp_data, exp_name in [(exp1_data, 'Exp1'), (exp2_data, 'Exp2'), (exp3_data, 'Exp3')]:
        for condition_name, condition_data in exp_data.items():
            tp_name, tp_hours = map_condition_to_timepoint(condition_name)

            if tp_hours is None:
                continue

            if tp_name not in pooled:
                pooled[tp_name] = {
                    'hours': tp_hours,
                    'c1': [],
                    'c2': [],
                    'is_treatment': 'treatment' in tp_name.lower(),
                    'is_recovery': 'recovery' in tp_name.lower()
                }

            pooled[tp_name]['c1'].extend(condition_data['c1'])
            pooled[tp_name]['c2'].extend(condition_data['c2'])

    return pooled

def plot_1_temporal_progression(pooled_data, results, output_dir):
    """Plot 1: Temporal progression with significance markers"""

    fig, ax = plt.subplots(figsize=(12, 6))

    # Separate controls and treatments
    controls = pooled_data['Controls']

    treatment_times = []
    treatment_means = []
    treatment_sems = []
    treatment_colors = []
    treatment_significant = []

    time_order = ['4h', '6h', '9h', '12h', '17h', '20h']

    for tp in time_order:
        if tp in pooled_data:
            data = pooled_data[tp]
            treatment_times.append(data['hours'])
            treatment_means.append(np.mean(data['c2']))
            treatment_sems.append(np.std(data['c2']) / np.sqrt(len(data['c2'])))

            # Check if significant
            result = next((r for r in results if r['time_point'] == tp), None)
            is_sig = result and result.get('significant_bonferroni', False)
            treatment_significant.append(is_sig)
            treatment_colors.append('#E74C3C' if is_sig else '#95A5A6')

    # Plot baseline (Controls)
    control_mean = np.mean(controls['c2'])
    control_sem = np.std(controls['c2']) / np.sqrt(len(controls['c2']))
    ax.axhline(control_mean, color='#3498DB', linestyle='--', linewidth=2,
               label='Controls (Baseline)', zorder=1)
    ax.axhspan(control_mean - control_sem, control_mean + control_sem,
               color='#3498DB', alpha=0.2, zorder=0)

    # Plot treatment progression
    ax.errorbar(treatment_times, treatment_means, yerr=treatment_sems,
                fmt='o-', linewidth=2, markersize=8, capsize=5,
                color='#E74C3C', label='Treatment', zorder=2)

    # Mark significant points
    for i, (time, mean, is_sig) in enumerate(zip(treatment_times, treatment_means, treatment_significant)):
        if is_sig:
            ax.plot(time, mean, 'o', markersize=12, markerfacecolor='#E74C3C',
                   markeredgecolor='#C0392B', markeredgewidth=3, zorder=3)
            ax.text(time, mean + 0.5, '***', ha='center', va='bottom',
                   fontsize=14, fontweight='bold', color='#C0392B')

    ax.set_xlabel('Treatment Duration (hours)', fontsize=12, fontweight='bold')
    ax.set_ylabel('C2 Scatter Index\n(Protein X Dispersion)', fontsize=12, fontweight='bold')
    ax.set_title('Drug Effect on Protein Scattering Over Time\n(Pooled Analysis: 3 Experiments, ~4,700 Cells)',
                fontsize=13, fontweight='bold', pad=15)
    ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    ax.set_xlim(-1, 21)

    # Add significance annotation
    ax.text(0.98, 0.02, '*** p < 0.008 (Bonferroni-corrected)',
           transform=ax.transAxes, ha='right', va='bottom',
           fontsize=9, style='italic', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_dir / 'plot1_temporal_progression.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: plot1_temporal_progression.png")
    plt.close()

def plot_2_effect_sizes(results, output_dir):
    """Plot 2: Effect sizes (Cohen's d) with significance"""

    fig, ax = plt.subplots(figsize=(10, 6))

    time_points = [r['time_point'] for r in results]
    cohens_d = [r['c2_cohens_d'] for r in results]
    is_significant = [r.get('significant_bonferroni', False) for r in results]

    colors = ['#E74C3C' if sig else '#95A5A6' for sig in is_significant]

    bars = ax.bar(range(len(time_points)), cohens_d, color=colors,
                  edgecolor='black', linewidth=1.5, alpha=0.8)

    # Add reference lines for effect size interpretation
    ax.axhline(0.2, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='Small effect (0.2)')
    ax.axhline(0.5, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='Medium effect (0.5)')
    ax.axhline(0.8, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Large effect (0.8)')

    # Mark significant bars
    for i, (bar, sig, d) in enumerate(zip(bars, is_significant, cohens_d)):
        if sig:
            ax.text(i, d + 0.03, '***', ha='center', va='bottom',
                   fontsize=14, fontweight='bold', color='#C0392B')

    ax.set_xticks(range(len(time_points)))
    ax.set_xticklabels(time_points, fontsize=11)
    ax.set_ylabel("Cohen's d (Effect Size)", fontsize=12, fontweight='bold')
    ax.set_xlabel('Treatment Duration', fontsize=12, fontweight='bold')
    ax.set_title('Effect Sizes: Drug Impact on Protein Scattering\n(Negative = decreased, Positive = increased)',
                fontsize=13, fontweight='bold', pad=15)
    ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True, fontsize=9)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5, axis='y')
    ax.set_ylim(-0.2, 0.6)

    plt.tight_layout()
    plt.savefig(output_dir / 'plot2_effect_sizes.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: plot2_effect_sizes.png")
    plt.close()

def plot_3_pvalue_heatmap(results, output_dir):
    """Plot 3: P-value significance heatmap"""

    fig, ax = plt.subplots(figsize=(10, 3))

    time_points = [r['time_point'] for r in results]
    p_values = [r['c2_p_value'] for r in results]

    # Create significance levels
    sig_levels = []
    for p in p_values:
        if p < 0.001:
            sig_levels.append(4)  # p < 0.001
        elif p < 0.008333:  # Bonferroni threshold
            sig_levels.append(3)  # p < 0.008
        elif p < 0.01:
            sig_levels.append(2)  # p < 0.01
        elif p < 0.05:
            sig_levels.append(1)  # p < 0.05
        else:
            sig_levels.append(0)  # ns

    # Create heatmap
    data = np.array([sig_levels])

    im = ax.imshow(data, cmap='RdYlGn', aspect='auto', vmin=0, vmax=4)

    ax.set_xticks(range(len(time_points)))
    ax.set_xticklabels(time_points, fontsize=11)
    ax.set_yticks([0])
    ax.set_yticklabels(['Statistical\nSignificance'], fontsize=11)

    # Add p-values as text
    for i, (tp, p, sig) in enumerate(zip(time_points, p_values, sig_levels)):
        color = 'white' if sig >= 3 else 'black'
        if p < 0.001:
            text = 'p<0.001\n***'
        else:
            text = f'p={p:.4f}\n{"***" if sig >= 3 else ("**" if sig == 2 else ("*" if sig == 1 else "ns"))}'
        ax.text(i, 0, text, ha='center', va='center', fontsize=9,
               fontweight='bold', color=color)

    ax.set_title('Statistical Significance Across Treatment Durations\n(Bonferroni-corrected threshold: p < 0.008)',
                fontsize=13, fontweight='bold', pad=15)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.15, aspect=30)
    cbar.set_ticks([0, 1, 2, 3, 4])
    cbar.set_ticklabels(['ns\n(p≥0.05)', 'p<0.05', 'p<0.01', 'p<0.008\n(Bonf.)', 'p<0.001'])
    cbar.ax.tick_params(labelsize=8)

    plt.tight_layout()
    plt.savefig(output_dir / 'plot3_pvalue_heatmap.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: plot3_pvalue_heatmap.png")
    plt.close()

def plot_4_scatter_distributions(pooled_data, output_dir):
    """Plot 4: Violin plots comparing distributions"""

    fig, ax = plt.subplots(figsize=(14, 6))

    # Prepare data
    time_order = ['Controls', '4h', '6h', '9h', '12h', '17h', '20h']
    plot_data = []
    plot_labels = []
    plot_colors = []

    for tp in time_order:
        if tp in pooled_data:
            plot_data.append(pooled_data[tp]['c2'])
            plot_labels.append(f"{tp}\n(n={len(pooled_data[tp]['c2'])})")
            plot_colors.append('#3498DB' if tp == 'Controls' else '#E74C3C')

    # Create violin plot
    parts = ax.violinplot(plot_data, positions=range(len(plot_data)),
                          showmeans=True, showmedians=True)

    # Color violins
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(plot_colors[i])
        pc.set_alpha(0.7)
        pc.set_edgecolor('black')
        pc.set_linewidth(1)

    # Style the other elements
    for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians', 'cmeans'):
        if partname in parts:
            parts[partname].set_edgecolor('black')
            parts[partname].set_linewidth(1.5)

    ax.set_xticks(range(len(plot_labels)))
    ax.set_xticklabels(plot_labels, fontsize=10)
    ax.set_ylabel('C2 Scatter Index (Protein X)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
    ax.set_title('Distribution of Protein Scattering Across Conditions\n(Violin plot shows full distribution + median + mean)',
                fontsize=13, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5, axis='y')

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#3498DB', alpha=0.7, edgecolor='black', label='Controls'),
        Patch(facecolor='#E74C3C', alpha=0.7, edgecolor='black', label='Treatment')
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True, fancybox=True, shadow=True)

    plt.tight_layout()
    plt.savefig(output_dir / 'plot4_scatter_distributions.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: plot4_scatter_distributions.png")
    plt.close()

def plot_5_sample_sizes(results, pooled_data, output_dir):
    """Plot 5: Sample sizes per condition"""

    fig, ax = plt.subplots(figsize=(10, 5))

    time_points = ['Controls'] + [r['time_point'] for r in results]
    sample_sizes = [len(pooled_data['Controls']['c2'])] + [r['n_treatment'] for r in results]

    colors = ['#3498DB'] + ['#E74C3C'] * len(results)

    bars = ax.bar(range(len(time_points)), sample_sizes, color=colors,
                  edgecolor='black', linewidth=1.5, alpha=0.8)

    # Add numbers on bars
    for i, (bar, size) in enumerate(zip(bars, sample_sizes)):
        ax.text(i, size + 10, str(size), ha='center', va='bottom',
               fontsize=10, fontweight='bold')

    ax.set_xticks(range(len(time_points)))
    ax.set_xticklabels(time_points, fontsize=11)
    ax.set_ylabel('Number of Cells Analyzed', fontsize=12, fontweight='bold')
    ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
    ax.set_title('Statistical Power: Sample Sizes Across Conditions\n(Pooled from 3 Independent Experiments)',
                fontsize=13, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5, axis='y')

    # Add total annotation
    total_cells = sum(sample_sizes)
    ax.text(0.98, 0.98, f'Total cells analyzed: {total_cells:,}',
           transform=ax.transAxes, ha='right', va='top',
           fontsize=11, fontweight='bold',
           bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    plt.tight_layout()
    plt.savefig(output_dir / 'plot5_sample_sizes.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: plot5_sample_sizes.png")
    plt.close()

def plot_6_percent_change_waterfall(results, output_dir):
    """Plot 6: Waterfall chart of percent changes"""

    fig, ax = plt.subplots(figsize=(10, 6))

    time_points = [r['time_point'] for r in results]
    percent_changes = [r['c2_percent_change'] for r in results]
    is_significant = [r.get('significant_bonferroni', False) for r in results]

    colors = ['#E74C3C' if sig else '#95A5A6' for sig in is_significant]

    # Create waterfall
    bars = ax.bar(range(len(time_points)), percent_changes, color=colors,
                  edgecolor='black', linewidth=1.5, alpha=0.8)

    # Add zero line
    ax.axhline(0, color='black', linestyle='-', linewidth=1.5)

    # Add values on bars
    for i, (bar, pct, sig) in enumerate(zip(bars, percent_changes, is_significant)):
        height = bar.get_height()
        ax.text(i, height + (1 if height > 0 else -1), f'{pct:+.1f}%',
               ha='center', va='bottom' if height > 0 else 'top',
               fontsize=10, fontweight='bold')
        if sig:
            ax.text(i, height + (2 if height > 0 else -2), '***',
                   ha='center', va='bottom' if height > 0 else 'top',
                   fontsize=14, fontweight='bold', color='#C0392B')

    ax.set_xticks(range(len(time_points)))
    ax.set_xticklabels(time_points, fontsize=11)
    ax.set_ylabel('Percent Change from Controls (%)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Treatment Duration', fontsize=12, fontweight='bold')
    ax.set_title('Relative Change in Protein Scattering vs Controls\n(Positive = Increased scattering)',
                fontsize=13, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5, axis='y')
    ax.set_ylim(-8, 20)

    plt.tight_layout()
    plt.savefig(output_dir / 'plot6_percent_change_waterfall.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: plot6_percent_change_waterfall.png")
    plt.close()

def main():
    print("\n" + "="*80)
    print("CREATING PUBLICATION-QUALITY PLOTS")
    print("="*80 + "\n")

    # Load results
    with open('results/pooled_analysis_results.json', 'r') as f:
        analysis_results = json.load(f)

    results = analysis_results['results']

    # Load raw data
    exp1_data, exp2_data, exp3_data = load_pooled_data()
    pooled_data = pool_by_timepoint(exp1_data, exp2_data, exp3_data)

    # Create output directory
    output_dir = Path('results/plots/png')
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\nGenerating plots...\n")

    # Generate all plots
    plot_1_temporal_progression(pooled_data, results, output_dir)
    plot_2_effect_sizes(results, output_dir)
    plot_3_pvalue_heatmap(results, output_dir)
    plot_4_scatter_distributions(pooled_data, output_dir)
    plot_5_sample_sizes(results, pooled_data, output_dir)
    plot_6_percent_change_waterfall(results, output_dir)

    print("\n" + "="*80)
    print(f"✅ All plots saved to: {output_dir}/")
    print("="*80)
    print("\nGenerated plots:")
    print("  1. plot1_temporal_progression.png - Main finding: scattering over time")
    print("  2. plot2_effect_sizes.png - Cohen's d effect sizes")
    print("  3. plot3_pvalue_heatmap.png - Statistical significance levels")
    print("  4. plot4_scatter_distributions.png - Full distributions (violin plots)")
    print("  5. plot5_sample_sizes.png - Sample sizes per condition")
    print("  6. plot6_percent_change_waterfall.png - Relative changes from baseline")
    print()

if __name__ == '__main__':
    main()
