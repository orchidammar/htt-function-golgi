#!/usr/bin/env python3
"""
Create plots showing both C1 (Structure Y) and C2 (Protein X) scatter individually
to demonstrate co-scattering over time, INCLUDING RECOVERY CONDITIONS.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set publication style
plt.style.use('seaborn-v0_8-paper')
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 10

OUTPUT_DIR = Path("results/plots/png/dual_channel")

def load_results():
    """Load pooled analysis results"""
    with open('results/pooled_analysis_results.json', 'r') as f:
        data = json.load(f)
    return data['results'], data['bonferroni_alpha']

def plot_dual_channel_temporal():
    """Plot C1 and C2 scatter over time INCLUDING RECOVERY"""
    results, bonferroni_alpha = load_results()

    # Treatment time points shown: 4h -> 6h -> 9h -> 20h (12h/17h excluded)
    treatment_labels = ['4h', '6h', '9h', '20h']
    recovery_labels = ['recovery_12h', 'recovery_17h', 'recovery_20h']

    # Categorical x positions: recovery sits AFTER the 20h treatment point
    treatment_hours = list(range(len(treatment_labels)))
    recovery_hours = list(range(len(treatment_labels),
                               len(treatment_labels) + len(recovery_labels)))
    tick_positions = treatment_hours + recovery_hours
    tick_labels = ['4h', '6h', '9h', '20h', 'Rec 12h', 'Rec 17h', 'Rec 20h']

    # Extract treatment data
    c1_treatment_means = []
    c1_treatment_stds = []
    c1_treatment_pvals = []
    c2_treatment_means = []
    c2_treatment_stds = []
    c2_treatment_pvals = []

    for tl in treatment_labels:
        r = next((r for r in results if r['time_point'] == tl), None)
        if r:
            c1_treatment_means.append(r['c1_treatment_mean'])
            c1_treatment_stds.append(r['c1_treatment_std'])
            c1_treatment_pvals.append(r['c1_p_value'])
            c2_treatment_means.append(r['c2_treatment_mean'])
            c2_treatment_stds.append(r['c2_treatment_std'])
            c2_treatment_pvals.append(r['c2_p_value'])

    # Extract recovery data
    c1_recovery_means = []
    c1_recovery_stds = []
    c1_recovery_pvals = []
    c2_recovery_means = []
    c2_recovery_stds = []
    c2_recovery_pvals = []

    for rl in recovery_labels:
        r = next((r for r in results if r['time_point'] == rl), None)
        if r:
            c1_recovery_means.append(r['c1_treatment_mean'])
            c1_recovery_stds.append(r['c1_treatment_std'])
            c1_recovery_pvals.append(r['c1_p_value'])
            c2_recovery_means.append(r['c2_treatment_mean'])
            c2_recovery_stds.append(r['c2_treatment_std'])
            c2_recovery_pvals.append(r['c2_p_value'])

    # Get baseline
    c1_baseline = results[0]['c1_control_mean']
    c2_baseline = results[0]['c2_control_mean']

    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot TREATMENT - C1 (Structure Y - RED) - solid line
    ax.errorbar(treatment_hours, c1_treatment_means, yerr=c1_treatment_stds,
                marker='o', markersize=8, linewidth=2.5, capsize=5,
                color='#d62728', label='C1: Structure Y (RED) - Treatment',
                alpha=0.8, linestyle='-')

    # Plot RECOVERY - C1 (Structure Y - RED) - dashed line
    ax.errorbar(recovery_hours, c1_recovery_means, yerr=c1_recovery_stds,
                marker='D', markersize=8, linewidth=2.5, capsize=5,
                color='#d62728', label='C1: Structure Y (RED) - RECOVERY',
                alpha=0.6, linestyle='--')

    # Plot TREATMENT - C2 (Protein X - GREEN) - solid line
    ax.errorbar(treatment_hours, c2_treatment_means, yerr=c2_treatment_stds,
                marker='s', markersize=8, linewidth=2.5, capsize=5,
                color='#2ca02c', label='C2: Protein X (GREEN) - Treatment',
                alpha=0.8, linestyle='-')

    # Plot RECOVERY - C2 (Protein X - GREEN) - dashed line
    ax.errorbar(recovery_hours, c2_recovery_means, yerr=c2_recovery_stds,
                marker='D', markersize=8, linewidth=2.5, capsize=5,
                color='#2ca02c', label='C2: Protein X (GREEN) - RECOVERY',
                alpha=0.6, linestyle='--')

    # Add baseline references
    ax.axhline(c1_baseline, color='#d62728', linestyle=':', alpha=0.4, linewidth=1.5)
    ax.axhline(c2_baseline, color='#2ca02c', linestyle=':', alpha=0.4, linewidth=1.5)

    # Mark significant treatment points with stars
    for i, (t, p1, p2) in enumerate(zip(treatment_hours, c1_treatment_pvals, c2_treatment_pvals)):
        if p1 < bonferroni_alpha:
            ax.plot(t, c1_treatment_means[i], marker='*', markersize=16, color='red', zorder=10)
        if p2 < bonferroni_alpha:
            ax.plot(t, c2_treatment_means[i], marker='*', markersize=16, color='darkgreen', zorder=10)

    # Mark significant recovery points with stars
    for i, (t, p1, p2) in enumerate(zip(recovery_hours, c1_recovery_pvals, c2_recovery_pvals)):
        if p1 < bonferroni_alpha:
            ax.plot(t, c1_recovery_means[i], marker='*', markersize=16, color='red', zorder=10)
        if p2 < bonferroni_alpha:
            ax.plot(t, c2_recovery_means[i], marker='*', markersize=16, color='darkgreen', zorder=10)

    ax.set_xlabel('Time point', fontsize=12, fontweight='bold')
    ax.set_ylabel('Scatter Index (a.u.)', fontsize=12, fontweight='bold')
    ax.set_title('Drug Effect and Recovery: Structure Y (RED) and Protein X (GREEN)\nSolid = Treatment | Dashed = RECOVERY',
                 fontsize=13, fontweight='bold', pad=15)
    ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True, fontsize=9)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)

    # Set x-axis ticks
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)

    # Add significance note
    ax.text(0.98, 0.02, f'* = p < {bonferroni_alpha:.4f} (Bonferroni-corrected)',
            transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'dual_channel_temporal_with_recovery.png', dpi=300, bbox_inches='tight')
    print(f"✅ Saved: {OUTPUT_DIR / 'dual_channel_temporal_with_recovery.png'}")
    plt.close()

def plot_percent_change_comparison():
    """Plot percent change from baseline INCLUDING RECOVERY"""
    results, bonferroni_alpha = load_results()

    # All time points including recovery
    all_labels = ['4h', '6h', '9h', '12h', '17h', '20h',
                  'Rec\n12h', 'Rec\n17h', 'Rec\n20h']
    all_keys = ['4h', '6h', '9h', '12h', '17h', '20h',
                'recovery_12h', 'recovery_17h', 'recovery_20h']

    c1_changes = []
    c2_changes = []
    c1_sig = []
    c2_sig = []

    for key in all_keys:
        r = next((r for r in results if r['time_point'] == key), None)
        if r:
            c1_changes.append(r['c1_percent_change'])
            c2_changes.append(r['c2_percent_change'])
            c1_sig.append(r.get('c1_significant_bonferroni', False))
            c2_sig.append(r.get('c2_significant_bonferroni', False))

    x = np.arange(len(all_labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 6))

    # Color-code treatment vs recovery
    c1_colors = ['#d62728']*6 + ['#ff9896']*3  # Darker for treatment, lighter for recovery
    c2_colors = ['#2ca02c']*6 + ['#98df8a']*3

    bars1 = ax.bar(x - width/2, c1_changes, width, label='C1: Structure Y (RED)',
                   color=c1_colors, alpha=0.7, edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width/2, c2_changes, width, label='C2: Protein X (GREEN)',
                   color=c2_colors, alpha=0.7, edgecolor='black', linewidth=1.5)

    # Mark significant bars
    for i, (b1, b2, s1, s2) in enumerate(zip(bars1, bars2, c1_sig, c2_sig)):
        if s1:
            height = b1.get_height()
            ax.text(b1.get_x() + b1.get_width()/2, height + 0.5 if height > 0 else height - 1.5,
                   '***', ha='center', va='bottom' if height > 0 else 'top',
                   fontsize=14, fontweight='bold', color='darkred')
        if s2:
            height = b2.get_height()
            ax.text(b2.get_x() + b2.get_width()/2, height + 0.5 if height > 0 else height - 1.5,
                   '***', ha='center', va='bottom' if height > 0 else 'top',
                   fontsize=14, fontweight='bold', color='darkgreen')

    # Add vertical line separating treatment and recovery
    ax.axvline(5.5, color='black', linestyle='--', alpha=0.5, linewidth=2)
    ax.text(5.5, ax.get_ylim()[1]*0.95, 'RECOVERY →', ha='left', va='top',
            fontsize=10, fontweight='bold', bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

    ax.axhline(0, color='black', linewidth=1.5, linestyle='-', alpha=0.5)
    ax.set_xlabel('Time Point', fontsize=12, fontweight='bold')
    ax.set_ylabel('Percent Change from Controls (%)', fontsize=12, fontweight='bold')
    ax.set_title('Scatter Change: Treatment and RECOVERY\nStructure Y (RED) vs Protein X (GREEN)',
                 fontsize=13, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(all_labels, fontsize=9)
    ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True)
    ax.grid(True, alpha=0.3, axis='y', linestyle=':', linewidth=0.5)

    # Add significance note
    ax.text(0.98, 0.02, f'*** = p < {bonferroni_alpha:.4f} (Bonferroni-corrected)',
            transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'dual_channel_percent_change_with_recovery.png', dpi=300, bbox_inches='tight')
    print(f"✅ Saved: {OUTPUT_DIR / 'dual_channel_percent_change_with_recovery.png'}")
    plt.close()

def plot_scatter_plot_c1_vs_c2():
    """Scatter plot showing C1 change vs C2 change INCLUDING RECOVERY"""
    results, bonferroni_alpha = load_results()

    # Separate treatment and recovery
    treatment_labels = ['4h', '6h', '9h', '12h', '17h', '20h']
    recovery_labels = ['recovery_12h', 'recovery_17h', 'recovery_20h']

    fig, ax = plt.subplots(figsize=(9, 9))

    # Plot treatment points
    for tl in treatment_labels:
        r = next((r for r in results if r['time_point'] == tl), None)
        if r:
            c1 = r['c1_percent_change']
            c2 = r['c2_percent_change']
            both_sig = r.get('c1_significant_bonferroni', False) and r.get('c2_significant_bonferroni', False)

            color = 'red' if both_sig else 'gray'
            marker = 'o'
            size = 200 if both_sig else 100
            ax.scatter(c1, c2, s=size, marker=marker, color=color, alpha=0.7,
                      edgecolor='black', linewidth=2, zorder=10, label=f'{tl}')
            ax.annotate(tl, (c1, c2), xytext=(5, 5), textcoords='offset points',
                       fontsize=10, fontweight='bold' if both_sig else 'normal')

    # Plot recovery points (different marker)
    for rl in recovery_labels:
        r = next((r for r in results if r['time_point'] == rl), None)
        if r:
            c1 = r['c1_percent_change']
            c2 = r['c2_percent_change']
            both_sig = r.get('c1_significant_bonferroni', False) and r.get('c2_significant_bonferroni', False)

            color = 'orange' if both_sig else 'lightgray'
            marker = 'D'  # Diamond for recovery
            size = 200 if both_sig else 100
            label = rl.replace('recovery_', 'Rec ')
            ax.scatter(c1, c2, s=size, marker=marker, color=color, alpha=0.7,
                      edgecolor='black', linewidth=2, zorder=10)
            ax.annotate(label, (c1, c2), xytext=(5, 5), textcoords='offset points',
                       fontsize=9, fontweight='bold' if both_sig else 'normal',
                       color='darkblue')

    # Add diagonal line (perfect correlation)
    all_c1 = [r['c1_percent_change'] for r in results]
    all_c2 = [r['c2_percent_change'] for r in results]
    min_val = min(min(all_c1), min(all_c2))
    max_val = max(max(all_c1), max(all_c2))
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.3, linewidth=2,
           label='Perfect correlation (1:1)')

    # Add quadrants
    ax.axhline(0, color='black', linewidth=1, alpha=0.5)
    ax.axvline(0, color='black', linewidth=1, alpha=0.5)

    ax.set_xlabel('C1 (Structure Y - RED) % Change', fontsize=12, fontweight='bold')
    ax.set_ylabel('C2 (Protein X - GREEN) % Change', fontsize=12, fontweight='bold')
    ax.set_title('Co-Scattering: Structure Y vs Protein X\nCircles = Treatment | Diamonds = RECOVERY',
                 fontsize=13, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    ax.set_aspect('equal', adjustable='box')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'dual_channel_correlation_with_recovery.png', dpi=300, bbox_inches='tight')
    print(f"✅ Saved: {OUTPUT_DIR / 'dual_channel_correlation_with_recovery.png'}")
    plt.close()

def plot_effect_size_comparison():
    """Compare Cohen's d effect sizes INCLUDING RECOVERY"""
    results, bonferroni_alpha = load_results()

    all_labels = ['4h', '6h', '9h', '12h', '17h', '20h',
                  'Rec\n12h', 'Rec\n17h', 'Rec\n20h']
    all_keys = ['4h', '6h', '9h', '12h', '17h', '20h',
                'recovery_12h', 'recovery_17h', 'recovery_20h']

    c1_d = []
    c2_d = []
    c1_sig = []
    c2_sig = []

    for key in all_keys:
        r = next((r for r in results if r['time_point'] == key), None)
        if r:
            c1_d.append(r['c1_cohens_d'])
            c2_d.append(r['c2_cohens_d'])
            c1_sig.append(r.get('c1_significant_bonferroni', False))
            c2_sig.append(r.get('c2_significant_bonferroni', False))

    x = np.arange(len(all_labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 6))

    # Color-code treatment vs recovery
    c1_colors = ['#d62728']*6 + ['#ff9896']*3
    c2_colors = ['#2ca02c']*6 + ['#98df8a']*3

    bars1 = ax.bar(x - width/2, c1_d, width, label='C1: Structure Y (RED)',
                   color=c1_colors, alpha=0.7, edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width/2, c2_d, width, label='C2: Protein X (GREEN)',
                   color=c2_colors, alpha=0.7, edgecolor='black', linewidth=1.5)

    # Mark significant bars
    for i, (b1, b2, s1, s2) in enumerate(zip(bars1, bars2, c1_sig, c2_sig)):
        if s1:
            height = b1.get_height()
            y_pos = height + 0.02 if height > 0 else height - 0.05
            ax.text(b1.get_x() + b1.get_width()/2, y_pos,
                   '***', ha='center', va='bottom' if height > 0 else 'top',
                   fontsize=14, fontweight='bold', color='darkred')
        if s2:
            height = b2.get_height()
            y_pos = height + 0.02 if height > 0 else height - 0.05
            ax.text(b2.get_x() + b2.get_width()/2, y_pos,
                   '***', ha='center', va='bottom' if height > 0 else 'top',
                   fontsize=14, fontweight='bold', color='darkgreen')

    # Add vertical line separating treatment and recovery
    ax.axvline(5.5, color='black', linestyle='--', alpha=0.5, linewidth=2)
    ax.text(5.5, ax.get_ylim()[1]*0.95, 'RECOVERY →', ha='left', va='top',
            fontsize=10, fontweight='bold', bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))

    # Add effect size reference lines
    ax.axhline(0.2, color='gray', linestyle=':', alpha=0.5, linewidth=1, label='Small effect (0.2)')
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1, label='Medium effect (0.5)')
    ax.axhline(0, color='black', linestyle='-', alpha=0.5, linewidth=1)

    ax.set_xlabel('Time Point', fontsize=12, fontweight='bold')
    ax.set_ylabel("Cohen's d (Effect Size)", fontsize=12, fontweight='bold')
    ax.set_title('Effect Size: Treatment and RECOVERY\nStructure Y (RED) vs Protein X (GREEN)',
                 fontsize=13, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(all_labels, fontsize=9)
    ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True, fontsize=9)
    ax.grid(True, alpha=0.3, axis='y', linestyle=':', linewidth=0.5)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'dual_channel_effect_sizes_with_recovery.png', dpi=300, bbox_inches='tight')
    print(f"✅ Saved: {OUTPUT_DIR / 'dual_channel_effect_sizes_with_recovery.png'}")
    plt.close()

def plot_reversibility_analysis():
    """Create dedicated reversibility plots showing recovery vs matched treatment"""
    results, bonferroni_alpha = load_results()

    # Match recovery to corresponding treatment duration
    pairs = [
        ('12h', 'recovery_12h', 12),
        ('17h', 'recovery_17h', 17),
        ('20h', 'recovery_20h', 20)
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for idx, (treatment_key, recovery_key, hours) in enumerate(pairs):
        ax = axes[idx]

        # Get data
        treatment = next((r for r in results if r['time_point'] == treatment_key), None)
        recovery = next((r for r in results if r['time_point'] == recovery_key), None)
        baseline_c1 = results[0]['c1_control_mean']
        baseline_c2 = results[0]['c2_control_mean']

        if treatment and recovery:
            # Data for plotting
            conditions = ['Controls', f'{hours}h\nTreatment', f'Recovery\n{hours}h']

            c1_values = [baseline_c1, treatment['c1_treatment_mean'], recovery['c1_treatment_mean']]
            c1_stds = [results[0]['c1_control_std'], treatment['c1_treatment_std'], recovery['c1_treatment_std']]

            c2_values = [baseline_c2, treatment['c2_treatment_mean'], recovery['c2_treatment_mean']]
            c2_stds = [results[0]['c2_control_std'], treatment['c2_treatment_std'], recovery['c2_treatment_std']]

            x = np.arange(len(conditions))
            width = 0.35

            # Plot bars
            bars1 = ax.bar(x - width/2, c1_values, width, yerr=c1_stds,
                          label='C1: Structure Y (RED)', color='#d62728', alpha=0.7,
                          capsize=5, edgecolor='black', linewidth=1.5)
            bars2 = ax.bar(x + width/2, c2_values, width, yerr=c2_stds,
                          label='C2: Protein X (GREEN)', color='#2ca02c', alpha=0.7,
                          capsize=5, edgecolor='black', linewidth=1.5)

            # Mark significance
            if treatment.get('c1_significant_bonferroni', False):
                ax.text(1 - width/2, treatment['c1_treatment_mean'] + treatment['c1_treatment_std'] + 0.5,
                       '***', ha='center', fontsize=12, color='darkred', fontweight='bold')
            if treatment.get('c2_significant_bonferroni', False):
                ax.text(1 + width/2, treatment['c2_treatment_mean'] + treatment['c2_treatment_std'] + 0.5,
                       '***', ha='center', fontsize=12, color='darkgreen', fontweight='bold')
            if recovery.get('c1_significant_bonferroni', False):
                y_pos = recovery['c1_treatment_mean'] + recovery['c1_treatment_std'] + 0.5
                ax.text(2 - width/2, y_pos,
                       '***', ha='center', fontsize=12, color='darkred', fontweight='bold')
            if recovery.get('c2_significant_bonferroni', False):
                y_pos = recovery['c2_treatment_mean'] + recovery['c2_treatment_std'] + 0.5
                ax.text(2 + width/2, y_pos,
                       '***', ha='center', fontsize=12, color='darkgreen', fontweight='bold')

            ax.set_ylabel('Scatter Index (a.u.)', fontsize=11, fontweight='bold')
            ax.set_title(f'{hours}h: Treatment → RECOVERY', fontsize=12, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels(conditions, fontsize=9)
            ax.grid(True, alpha=0.3, axis='y', linestyle=':', linewidth=0.5)

            if idx == 0:
                ax.legend(loc='upper left', fontsize=9)

    fig.suptitle('Reversibility Analysis: Does Scatter Return to Baseline After Drug Washout?',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'reversibility_analysis.png', dpi=300, bbox_inches='tight')
    print(f"✅ Saved: {OUTPUT_DIR / 'reversibility_analysis.png'}")
    plt.close()

def main():
    """Generate all dual-channel plots WITH RECOVERY"""
    print("\n" + "="*80)
    print("Creating dual-channel plots WITH RECOVERY DATA...")
    print("="*80)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n1. Temporal progression (both channels WITH RECOVERY)...")
    plot_dual_channel_temporal()

    print("2. Percent change comparison (bar chart WITH RECOVERY)...")
    plot_percent_change_comparison()

    print("3. C1 vs C2 correlation scatter plot (WITH RECOVERY)...")
    plot_scatter_plot_c1_vs_c2()

    print("4. Effect size comparison (WITH RECOVERY)...")
    plot_effect_size_comparison()

    print("5. Dedicated reversibility analysis...")
    plot_reversibility_analysis()

    print("\n" + "="*80)
    print(f"✅ All dual-channel plots WITH RECOVERY saved to: {OUTPUT_DIR}/")
    print("="*80)
    print("\nGenerated plots:")
    print("  1. dual_channel_temporal_with_recovery.png - Both channels over time + recovery")
    print("  2. dual_channel_percent_change_with_recovery.png - % change + recovery")
    print("  3. dual_channel_correlation_with_recovery.png - C1 vs C2 scatter + recovery")
    print("  4. dual_channel_effect_sizes_with_recovery.png - Cohen's d + recovery")
    print("  5. reversibility_analysis.png - Dedicated recovery comparison")
    print()

if __name__ == '__main__':
    main()
