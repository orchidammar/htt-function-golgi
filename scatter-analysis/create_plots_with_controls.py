#!/usr/bin/env python3
"""
Create alternative plots that include Controls as explicit data points
(not just baseline reference)
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set publication style
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("Set2")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 9

# Output directory for plots with controls
OUTPUT_DIR = Path("results/plots/png/plots_with_controls")

def load_pooled_data():
    """Load scatter data from all three experiments"""
    def load_from_experiment(output_dir):
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
                    'c2': []
                }

            pooled[tp_name]['c1'].extend(condition_data['c1'])
            pooled[tp_name]['c2'].extend(condition_data['c2'])

    return pooled

def plot_alt1_temporal_with_controls(pooled_data, results, output_dir):
    """Alternative Plot 1: Temporal progression WITH controls at time=0"""

    fig, ax = plt.subplots(figsize=(12, 6))

    # Include controls at time 0
    all_conditions = ['Controls'] + ['4h', '6h', '9h', '12h', '17h', '20h']

    times = []
    means = []
    sems = []
    colors = []
    markers = []
    significant = []

    for i, cond in enumerate(all_conditions):
        if cond in pooled_data:
            data = pooled_data[cond]
            times.append(data['hours'])
            means.append(np.mean(data['c2']))
            sems.append(np.std(data['c2']) / np.sqrt(len(data['c2'])))

            if cond == 'Controls':
                colors.append('#3498DB')
                markers.append('s')  # Square for controls
                significant.append(False)
            else:
                # Check if significant
                result = next((r for r in results if r['time_point'] == cond), None)
                is_sig = result and result.get('significant_bonferroni', False)
                significant.append(is_sig)
                colors.append('#E74C3C' if is_sig else '#95A5A6')
                markers.append('o')

    # Plot line
    ax.plot(times, means, 'o-', linewidth=2.5, color='#34495E',
            markersize=0, alpha=0.3, zorder=1)

    # Plot points with different colors
    for i, (time, mean, sem, color, marker, sig) in enumerate(zip(times, means, sems, colors, markers, significant)):
        # Error bar
        ax.errorbar(time, mean, yerr=sem, fmt=marker, markersize=10,
                   color=color, markeredgecolor='black', markeredgewidth=2,
                   capsize=5, capthick=2, elinewidth=2, zorder=2)

        # Mark significant points with stars
        if sig:
            ax.text(time, mean + sem + 0.5, '***', ha='center', va='bottom',
                   fontsize=14, fontweight='bold', color='#C0392B')

        # Label controls
        if i == 0:
            ax.text(time, mean - sem - 0.8, 'Baseline\n(Controls)', ha='center', va='top',
                   fontsize=9, style='italic', color='#2C3E50',
                   bbox=dict(boxstyle='round,pad=0.4', facecolor='#ECF0F1', alpha=0.8))

    ax.set_xlabel('Time (hours)', fontsize=12, fontweight='bold')
    ax.set_ylabel('C2 Scatter Index\n(Protein X Dispersion)', fontsize=12, fontweight='bold')
    ax.set_title('Drug Effect on Protein Scattering: Time Course with Baseline\n(Pooled Analysis: 3 Experiments, 4,449 Cells)',
                fontsize=13, fontweight='bold', pad=15)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='s', color='w', markerfacecolor='#3498DB',
               markeredgecolor='black', markersize=10, label='Controls (Baseline)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#E74C3C',
               markeredgecolor='black', markersize=10, label='Treatment (Significant)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#95A5A6',
               markeredgecolor='black', markersize=10, label='Treatment (Not Sig.)')
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True, fancybox=True, shadow=True)

    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    ax.set_xlim(-1, 21)
    ax.set_ylim(10, 15)

    # Add significance annotation
    ax.text(0.98, 0.02, '*** p < 0.008 (Bonferroni-corrected)',
           transform=ax.transAxes, ha='right', va='bottom',
           fontsize=9, style='italic', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_dir / 'plot_alt1_temporal_with_controls.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: plot_alt1_temporal_with_controls.png")
    plt.close()

def plot_alt2_bar_chart_all_conditions(pooled_data, results, output_dir):
    """Alternative Plot 2: Bar chart including controls"""

    fig, ax = plt.subplots(figsize=(12, 6))

    all_conditions = ['Controls', '4h', '6h', '9h', '12h', '17h', '20h']

    means = []
    sems = []
    colors = []
    labels = []
    n_values = []

    for cond in all_conditions:
        if cond in pooled_data:
            data = pooled_data[cond]
            means.append(np.mean(data['c2']))
            sems.append(np.std(data['c2']) / np.sqrt(len(data['c2'])))
            n_values.append(len(data['c2']))

            if cond == 'Controls':
                colors.append('#3498DB')
                labels.append(f'Controls\n(n={len(data["c2"])})')
            else:
                result = next((r for r in results if r['time_point'] == cond), None)
                is_sig = result and result.get('significant_bonferroni', False)
                colors.append('#E74C3C' if is_sig else '#95A5A6')
                labels.append(f'{cond}\n(n={len(data["c2"])})')

    x_pos = np.arange(len(labels))
    bars = ax.bar(x_pos, means, yerr=sems, color=colors,
                  edgecolor='black', linewidth=1.5, alpha=0.8,
                  capsize=5, error_kw={'elinewidth': 2, 'capthick': 2})

    # Add values on top of bars
    for i, (bar, mean, sem) in enumerate(zip(bars, means, sems)):
        ax.text(i, mean + sem + 0.2, f'{mean:.2f}', ha='center', va='bottom',
               fontsize=9, fontweight='bold')

        # Add significance stars
        if i > 0:  # Skip controls
            result = results[i-1] if i-1 < len(results) else None
            if result and result.get('significant_bonferroni', False):
                ax.text(i, mean + sem + 0.6, '***', ha='center', va='bottom',
                       fontsize=14, fontweight='bold', color='#C0392B')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel('C2 Scatter Index\n(Protein X Dispersion)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
    ax.set_title('Protein Scattering Across All Conditions\n(Mean ± SEM, *** = p<0.008)',
                fontsize=13, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5, axis='y')
    ax.set_ylim(10, 15)

    # Add baseline reference line
    baseline = means[0]
    ax.axhline(baseline, color='#3498DB', linestyle='--', linewidth=2, alpha=0.5,
              label=f'Baseline (Controls: {baseline:.2f})')
    ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True)

    plt.tight_layout()
    plt.savefig(output_dir / 'plot_alt2_bar_chart_all_conditions.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: plot_alt2_bar_chart_all_conditions.png")
    plt.close()

def plot_alt3_side_by_side_comparison(pooled_data, results, output_dir):
    """Alternative Plot 3: Side-by-side comparison showing absolute values"""

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    all_conditions = ['Controls', '4h', '6h', '9h', '12h', '17h', '20h']

    # LEFT PLOT: Absolute values
    means = []
    sems = []
    colors = []

    for cond in all_conditions:
        if cond in pooled_data:
            data = pooled_data[cond]
            means.append(np.mean(data['c2']))
            sems.append(np.std(data['c2']) / np.sqrt(len(data['c2'])))

            if cond == 'Controls':
                colors.append('#3498DB')
            else:
                result = next((r for r in results if r['time_point'] == cond), None)
                is_sig = result and result.get('significant_bonferroni', False)
                colors.append('#E74C3C' if is_sig else '#95A5A6')

    x_pos = np.arange(len(all_conditions))
    bars1 = ax1.bar(x_pos, means, yerr=sems, color=colors,
                    edgecolor='black', linewidth=1.5, alpha=0.8,
                    capsize=5, error_kw={'elinewidth': 2})

    for i, (bar, mean, sem) in enumerate(zip(bars1, means, sems)):
        if i > 0:
            result = results[i-1] if i-1 < len(results) else None
            if result and result.get('significant_bonferroni', False):
                ax1.text(i, mean + sem + 0.2, '***', ha='center', va='bottom',
                        fontsize=12, fontweight='bold', color='#C0392B')

    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(all_conditions, fontsize=10)
    ax1.set_ylabel('C2 Scatter Index', fontsize=11, fontweight='bold')
    ax1.set_xlabel('Condition', fontsize=11, fontweight='bold')
    ax1.set_title('A) Absolute Scatter Values\n(Mean ± SEM)', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle=':', linewidth=0.5, axis='y')
    ax1.axhline(means[0], color='#3498DB', linestyle='--', linewidth=1.5, alpha=0.5)
    ax1.set_ylim(10, 15)

    # RIGHT PLOT: Relative change from controls
    percent_changes = []
    for i, mean in enumerate(means):
        if i == 0:
            percent_changes.append(0)  # Controls = 0% change
        else:
            pct_change = ((mean - means[0]) / means[0]) * 100
            percent_changes.append(pct_change)

    colors_right = ['#3498DB'] + colors[1:]  # Controls blue, rest same as left
    bars2 = ax2.bar(x_pos, percent_changes, color=colors_right,
                    edgecolor='black', linewidth=1.5, alpha=0.8)

    ax2.axhline(0, color='black', linestyle='-', linewidth=2)

    for i, (bar, pct) in enumerate(zip(bars2, percent_changes)):
        if i > 0:
            result = results[i-1] if i-1 < len(results) else None
            if result and result.get('significant_bonferroni', False):
                ax2.text(i, pct + (1 if pct > 0 else -1), '***', ha='center',
                        va='bottom' if pct > 0 else 'top',
                        fontsize=12, fontweight='bold', color='#C0392B')
            ax2.text(i, pct/2, f'{pct:+.1f}%', ha='center', va='center',
                    fontsize=9, fontweight='bold', color='white' if abs(pct) > 5 else 'black')

    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(all_conditions, fontsize=10)
    ax2.set_ylabel('Percent Change from Controls (%)', fontsize=11, fontweight='bold')
    ax2.set_xlabel('Condition', fontsize=11, fontweight='bold')
    ax2.set_title('B) Relative Change from Baseline\n(Controls = 0%)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle=':', linewidth=0.5, axis='y')
    ax2.set_ylim(-8, 20)

    plt.suptitle('Protein Scattering: Absolute vs Relative Comparison\n(*** p < 0.008, Bonferroni-corrected)',
                fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout()
    plt.savefig(output_dir / 'plot_alt3_side_by_side_comparison.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: plot_alt3_side_by_side_comparison.png")
    plt.close()

def plot_alt4_grouped_bar_with_stats(pooled_data, results, output_dir):
    """Alternative Plot 4: Grouped bars with statistical annotations"""

    fig, ax = plt.subplots(figsize=(14, 7))

    all_conditions = ['Controls', '4h', '6h', '9h', '12h', '17h', '20h']

    x = np.arange(len(all_conditions))
    width = 0.6

    means = []
    sems = []
    p_values_text = []

    for i, cond in enumerate(all_conditions):
        if cond in pooled_data:
            data = pooled_data[cond]
            means.append(np.mean(data['c2']))
            sems.append(np.std(data['c2']) / np.sqrt(len(data['c2'])))

            if cond == 'Controls':
                p_values_text.append('Baseline')
            else:
                result = next((r for r in results if r['time_point'] == cond), None)
                if result:
                    p = result['c2_p_value']
                    if p < 0.001:
                        p_values_text.append('p<0.001***')
                    elif p < 0.008333:
                        p_values_text.append(f'p={p:.4f}***')
                    elif p < 0.05:
                        p_values_text.append(f'p={p:.3f}*')
                    else:
                        p_values_text.append(f'p={p:.3f}')
                else:
                    p_values_text.append('N/A')

    # Color code
    colors = []
    for i, cond in enumerate(all_conditions):
        if cond == 'Controls':
            colors.append('#3498DB')
        else:
            result = next((r for r in results if r['time_point'] == cond), None)
            is_sig = result and result.get('significant_bonferroni', False)
            colors.append('#E74C3C' if is_sig else '#95A5A6')

    bars = ax.bar(x, means, width, yerr=sems, color=colors,
                  edgecolor='black', linewidth=2, alpha=0.85,
                  capsize=6, error_kw={'elinewidth': 2.5, 'capthick': 2})

    # Add mean values on bars
    for i, (bar, mean, sem, p_text) in enumerate(zip(bars, means, sems, p_values_text)):
        # Mean value
        ax.text(i, mean + sem + 0.3, f'{mean:.2f}', ha='center', va='bottom',
               fontsize=10, fontweight='bold')

        # P-value below bar
        ax.text(i, 9.5, p_text, ha='center', va='top',
               fontsize=8, style='italic',
               color='#2C3E50' if i == 0 else ('#C0392B' if '***' in p_text else '#7F8C8D'))

        # Sample size
        n = len(pooled_data[all_conditions[i]]['c2'])
        ax.text(i, 9.2, f'n={n}', ha='center', va='top',
               fontsize=7, color='#34495E')

    ax.set_xticks(x)
    ax.set_xticklabels(all_conditions, fontsize=11, fontweight='bold')
    ax.set_ylabel('C2 Scatter Index (Protein X)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
    ax.set_title('Comprehensive Comparison: All Conditions with Statistics\n(Bars = Mean ± SEM, *** = Bonferroni-corrected p<0.008)',
                fontsize=13, fontweight='bold', pad=20)

    # Baseline reference
    baseline = means[0]
    ax.axhline(baseline, color='#3498DB', linestyle='--', linewidth=2, alpha=0.4,
              label=f'Controls Baseline: {baseline:.2f}')

    ax.grid(True, alpha=0.2, linestyle=':', linewidth=0.5, axis='y')
    ax.set_ylim(9, 15)
    ax.legend(loc='upper left', fontsize=10, frameon=True, fancybox=True, shadow=True)

    # Add significance zone
    ax.axhspan(baseline - sems[0], baseline + sems[0], color='#3498DB', alpha=0.1, zorder=0)

    plt.tight_layout()
    plt.savefig(output_dir / 'plot_alt4_grouped_bar_with_stats.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved: plot_alt4_grouped_bar_with_stats.png")
    plt.close()

def main():
    print("\n" + "="*80)
    print("CREATING ALTERNATIVE PLOTS (WITH CONTROLS VISIBLE)")
    print("="*80 + "\n")

    # Load results
    with open('results/pooled_analysis_results.json', 'r') as f:
        analysis_results = json.load(f)

    results = analysis_results['results']

    # Load raw data
    exp1_data, exp2_data, exp3_data = load_pooled_data()
    pooled_data = pool_by_timepoint(exp1_data, exp2_data, exp3_data)

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\nGenerating alternative plots...\n")

    # Generate alternative plots
    plot_alt1_temporal_with_controls(pooled_data, results, OUTPUT_DIR)
    plot_alt2_bar_chart_all_conditions(pooled_data, results, OUTPUT_DIR)
    plot_alt3_side_by_side_comparison(pooled_data, results, OUTPUT_DIR)
    plot_alt4_grouped_bar_with_stats(pooled_data, results, OUTPUT_DIR)

    print("\n" + "="*80)
    print(f"✅ All alternative plots saved to: {OUTPUT_DIR}/")
    print("="*80)
    print("\nGenerated alternative plots:")
    print("  1. plot_alt1_temporal_with_controls.png - Time course WITH controls at t=0")
    print("  2. plot_alt2_bar_chart_all_conditions.png - Bar chart INCLUDING controls")
    print("  3. plot_alt3_side_by_side_comparison.png - Absolute vs Relative (side-by-side)")
    print("  4. plot_alt4_grouped_bar_with_stats.png - All conditions with p-values annotated")
    print()

if __name__ == '__main__':
    main()
