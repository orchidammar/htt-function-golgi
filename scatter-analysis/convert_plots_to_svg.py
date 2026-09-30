#!/usr/bin/env python3
"""
Convert all plots to SVG format for vector editing.
SVG files can be edited in Adobe Illustrator, Inkscape, Figma, and even PowerPoint.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import seaborn as sns
from pathlib import Path
import os

# Set publication style with vector-friendly settings
plt.style.use('seaborn-v0_8-paper')
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 100
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['svg.fonttype'] = 'none'  # Keep text as text (editable)

# Output directory
SVG_DIR = Path("results/plots/svg")
SVG_DUAL_DIR = SVG_DIR / "dual_channel"

def load_results():
    """Load pooled analysis results"""
    with open('results/pooled_analysis_results.json', 'r') as f:
        data = json.load(f)
    return data['results'], data['bonferroni_alpha']

# ============================================================================
# DUAL CHANNEL PLOTS
# ============================================================================

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

    c1_treatment_means, c1_treatment_stds, c1_treatment_pvals = [], [], []
    c2_treatment_means, c2_treatment_stds, c2_treatment_pvals = [], [], []

    for tl in treatment_labels:
        r = next((r for r in results if r['time_point'] == tl), None)
        if r:
            c1_treatment_means.append(r['c1_treatment_mean'])
            c1_treatment_stds.append(r['c1_treatment_std'])
            c1_treatment_pvals.append(r['c1_p_value'])
            c2_treatment_means.append(r['c2_treatment_mean'])
            c2_treatment_stds.append(r['c2_treatment_std'])
            c2_treatment_pvals.append(r['c2_p_value'])

    c1_recovery_means, c1_recovery_stds, c1_recovery_pvals = [], [], []
    c2_recovery_means, c2_recovery_stds, c2_recovery_pvals = [], [], []

    for rl in recovery_labels:
        r = next((r for r in results if r['time_point'] == rl), None)
        if r:
            c1_recovery_means.append(r['c1_treatment_mean'])
            c1_recovery_stds.append(r['c1_treatment_std'])
            c1_recovery_pvals.append(r['c1_p_value'])
            c2_recovery_means.append(r['c2_treatment_mean'])
            c2_recovery_stds.append(r['c2_treatment_std'])
            c2_recovery_pvals.append(r['c2_p_value'])

    c1_baseline = results[0]['c1_control_mean']
    c2_baseline = results[0]['c2_control_mean']

    fig, ax = plt.subplots(figsize=(12, 6))

    ax.errorbar(treatment_hours, c1_treatment_means, yerr=c1_treatment_stds,
                marker='o', markersize=8, linewidth=2.5, capsize=5,
                color='#d62728', label='C1: Structure Y (RED) - Treatment',
                alpha=0.8, linestyle='-')
    ax.errorbar(recovery_hours, c1_recovery_means, yerr=c1_recovery_stds,
                marker='D', markersize=8, linewidth=2.5, capsize=5,
                color='#d62728', label='C1: Structure Y (RED) - RECOVERY',
                alpha=0.6, linestyle='--')
    ax.errorbar(treatment_hours, c2_treatment_means, yerr=c2_treatment_stds,
                marker='s', markersize=8, linewidth=2.5, capsize=5,
                color='#2ca02c', label='C2: Protein X (GREEN) - Treatment',
                alpha=0.8, linestyle='-')
    ax.errorbar(recovery_hours, c2_recovery_means, yerr=c2_recovery_stds,
                marker='D', markersize=8, linewidth=2.5, capsize=5,
                color='#2ca02c', label='C2: Protein X (GREEN) - RECOVERY',
                alpha=0.6, linestyle='--')

    ax.axhline(c1_baseline, color='#d62728', linestyle=':', alpha=0.4, linewidth=1.5)
    ax.axhline(c2_baseline, color='#2ca02c', linestyle=':', alpha=0.4, linewidth=1.5)

    for i, (t, p1, p2) in enumerate(zip(treatment_hours, c1_treatment_pvals, c2_treatment_pvals)):
        if p1 < bonferroni_alpha:
            ax.plot(t, c1_treatment_means[i], marker='*', markersize=16, color='red', zorder=10)
        if p2 < bonferroni_alpha:
            ax.plot(t, c2_treatment_means[i], marker='*', markersize=16, color='darkgreen', zorder=10)

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
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.text(0.98, 0.02, f'* = p < {bonferroni_alpha:.4f} (Bonferroni-corrected)',
            transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()
    plt.savefig(SVG_DUAL_DIR / 'dual_channel_temporal_with_recovery.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ dual_channel_temporal_with_recovery.svg")
    plt.close()


def plot_percent_change_comparison():
    """Plot percent change from baseline INCLUDING RECOVERY"""
    results, bonferroni_alpha = load_results()

    all_labels = ['4h', '6h', '9h', '12h', '17h', '20h', 'Rec\n12h', 'Rec\n17h', 'Rec\n20h']
    all_keys = ['4h', '6h', '9h', '12h', '17h', '20h', 'recovery_12h', 'recovery_17h', 'recovery_20h']

    c1_changes, c2_changes, c1_sig, c2_sig = [], [], [], []

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

    c1_colors = ['#d62728']*6 + ['#ff9896']*3
    c2_colors = ['#2ca02c']*6 + ['#98df8a']*3

    bars1 = ax.bar(x - width/2, c1_changes, width, label='C1: Structure Y (RED)',
                   color=c1_colors, alpha=0.7, edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width/2, c2_changes, width, label='C2: Protein X (GREEN)',
                   color=c2_colors, alpha=0.7, edgecolor='black', linewidth=1.5)

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
    ax.text(0.98, 0.02, f'*** = p < {bonferroni_alpha:.4f} (Bonferroni-corrected)',
            transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()
    plt.savefig(SVG_DUAL_DIR / 'dual_channel_percent_change_with_recovery.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ dual_channel_percent_change_with_recovery.svg")
    plt.close()


def plot_scatter_plot_c1_vs_c2():
    """Scatter plot showing C1 change vs C2 change INCLUDING RECOVERY"""
    results, bonferroni_alpha = load_results()

    treatment_labels = ['4h', '6h', '9h', '12h', '17h', '20h']
    recovery_labels = ['recovery_12h', 'recovery_17h', 'recovery_20h']

    fig, ax = plt.subplots(figsize=(9, 9))

    for tl in treatment_labels:
        r = next((r for r in results if r['time_point'] == tl), None)
        if r:
            c1 = r['c1_percent_change']
            c2 = r['c2_percent_change']
            both_sig = r.get('c1_significant_bonferroni', False) and r.get('c2_significant_bonferroni', False)
            color = 'red' if both_sig else 'gray'
            size = 200 if both_sig else 100
            ax.scatter(c1, c2, s=size, marker='o', color=color, alpha=0.7,
                      edgecolor='black', linewidth=2, zorder=10)
            ax.annotate(tl, (c1, c2), xytext=(5, 5), textcoords='offset points',
                       fontsize=10, fontweight='bold' if both_sig else 'normal')

    for rl in recovery_labels:
        r = next((r for r in results if r['time_point'] == rl), None)
        if r:
            c1 = r['c1_percent_change']
            c2 = r['c2_percent_change']
            both_sig = r.get('c1_significant_bonferroni', False) and r.get('c2_significant_bonferroni', False)
            color = 'orange' if both_sig else 'lightgray'
            size = 200 if both_sig else 100
            label = rl.replace('recovery_', 'Rec ')
            ax.scatter(c1, c2, s=size, marker='D', color=color, alpha=0.7,
                      edgecolor='black', linewidth=2, zorder=10)
            ax.annotate(label, (c1, c2), xytext=(5, 5), textcoords='offset points',
                       fontsize=9, fontweight='bold' if both_sig else 'normal', color='darkblue')

    all_c1 = [r['c1_percent_change'] for r in results]
    all_c2 = [r['c2_percent_change'] for r in results]
    min_val = min(min(all_c1), min(all_c2))
    max_val = max(max(all_c1), max(all_c2))
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.3, linewidth=2,
           label='Perfect correlation (1:1)')
    ax.axhline(0, color='black', linewidth=1, alpha=0.5)
    ax.axvline(0, color='black', linewidth=1, alpha=0.5)
    ax.set_xlabel('C1 (Structure Y - RED) % Change', fontsize=12, fontweight='bold')
    ax.set_ylabel('C2 (Protein X - GREEN) % Change', fontsize=12, fontweight='bold')
    ax.set_title('Co-Scattering: Structure Y vs Protein X\nCircles = Treatment | Diamonds = RECOVERY',
                 fontsize=13, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    ax.set_aspect('equal', adjustable='box')

    plt.tight_layout()
    plt.savefig(SVG_DUAL_DIR / 'dual_channel_correlation_with_recovery.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ dual_channel_correlation_with_recovery.svg")
    plt.close()


def plot_effect_size_comparison():
    """Compare Cohen's d effect sizes INCLUDING RECOVERY"""
    results, bonferroni_alpha = load_results()

    all_labels = ['4h', '6h', '9h', '12h', '17h', '20h', 'Rec\n12h', 'Rec\n17h', 'Rec\n20h']
    all_keys = ['4h', '6h', '9h', '12h', '17h', '20h', 'recovery_12h', 'recovery_17h', 'recovery_20h']

    c1_d, c2_d, c1_sig, c2_sig = [], [], [], []

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

    c1_colors = ['#d62728']*6 + ['#ff9896']*3
    c2_colors = ['#2ca02c']*6 + ['#98df8a']*3

    bars1 = ax.bar(x - width/2, c1_d, width, label='C1: Structure Y (RED)',
                   color=c1_colors, alpha=0.7, edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width/2, c2_d, width, label='C2: Protein X (GREEN)',
                   color=c2_colors, alpha=0.7, edgecolor='black', linewidth=1.5)

    for i, (b1, b2, s1, s2) in enumerate(zip(bars1, bars2, c1_sig, c2_sig)):
        if s1:
            height = b1.get_height()
            y_pos = height + 0.02 if height > 0 else height - 0.05
            ax.text(b1.get_x() + b1.get_width()/2, y_pos, '***', ha='center',
                   va='bottom' if height > 0 else 'top', fontsize=14, fontweight='bold', color='darkred')
        if s2:
            height = b2.get_height()
            y_pos = height + 0.02 if height > 0 else height - 0.05
            ax.text(b2.get_x() + b2.get_width()/2, y_pos, '***', ha='center',
                   va='bottom' if height > 0 else 'top', fontsize=14, fontweight='bold', color='darkgreen')

    ax.axvline(5.5, color='black', linestyle='--', alpha=0.5, linewidth=2)
    ax.text(5.5, ax.get_ylim()[1]*0.95, 'RECOVERY →', ha='left', va='top',
            fontsize=10, fontweight='bold', bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))
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
    plt.savefig(SVG_DUAL_DIR / 'dual_channel_effect_sizes_with_recovery.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ dual_channel_effect_sizes_with_recovery.svg")
    plt.close()


def plot_reversibility_analysis():
    """Create dedicated reversibility plots"""
    results, bonferroni_alpha = load_results()

    pairs = [('12h', 'recovery_12h', 12), ('17h', 'recovery_17h', 17), ('20h', 'recovery_20h', 20)]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for idx, (treatment_key, recovery_key, hours) in enumerate(pairs):
        ax = axes[idx]
        treatment = next((r for r in results if r['time_point'] == treatment_key), None)
        recovery = next((r for r in results if r['time_point'] == recovery_key), None)
        baseline_c1 = results[0]['c1_control_mean']
        baseline_c2 = results[0]['c2_control_mean']

        if treatment and recovery:
            conditions = ['Controls', f'{hours}h\nTreatment', f'Recovery\n{hours}h']
            c1_values = [baseline_c1, treatment['c1_treatment_mean'], recovery['c1_treatment_mean']]
            c1_stds = [results[0]['c1_control_std'], treatment['c1_treatment_std'], recovery['c1_treatment_std']]
            c2_values = [baseline_c2, treatment['c2_treatment_mean'], recovery['c2_treatment_mean']]
            c2_stds = [results[0]['c2_control_std'], treatment['c2_treatment_std'], recovery['c2_treatment_std']]

            x = np.arange(len(conditions))
            width = 0.35

            ax.bar(x - width/2, c1_values, width, yerr=c1_stds, label='C1: Structure Y (RED)',
                  color='#d62728', alpha=0.7, capsize=5, edgecolor='black', linewidth=1.5)
            ax.bar(x + width/2, c2_values, width, yerr=c2_stds, label='C2: Protein X (GREEN)',
                  color='#2ca02c', alpha=0.7, capsize=5, edgecolor='black', linewidth=1.5)

            if treatment.get('c1_significant_bonferroni', False):
                ax.text(1 - width/2, treatment['c1_treatment_mean'] + treatment['c1_treatment_std'] + 0.5,
                       '***', ha='center', fontsize=12, color='darkred', fontweight='bold')
            if treatment.get('c2_significant_bonferroni', False):
                ax.text(1 + width/2, treatment['c2_treatment_mean'] + treatment['c2_treatment_std'] + 0.5,
                       '***', ha='center', fontsize=12, color='darkgreen', fontweight='bold')
            if recovery.get('c1_significant_bonferroni', False):
                ax.text(2 - width/2, recovery['c1_treatment_mean'] + recovery['c1_treatment_std'] + 0.5,
                       '***', ha='center', fontsize=12, color='darkred', fontweight='bold')
            if recovery.get('c2_significant_bonferroni', False):
                ax.text(2 + width/2, recovery['c2_treatment_mean'] + recovery['c2_treatment_std'] + 0.5,
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
    plt.savefig(SVG_DUAL_DIR / 'reversibility_analysis.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ reversibility_analysis.svg")
    plt.close()


# ============================================================================
# SINGLE CHANNEL PLOTS (from create_pooled_plots.py)
# ============================================================================

def plot_temporal_progression():
    """C2 scatter over treatment time"""
    results, bonferroni_alpha = load_results()

    treatment_labels = ['4h', '6h', '9h', '12h', '17h', '20h']
    treatment_hours = [4, 6, 9, 12, 17, 20]

    means, stds, pvals = [], [], []
    for tl in treatment_labels:
        r = next((r for r in results if r['time_point'] == tl), None)
        if r:
            means.append(r['c2_treatment_mean'])
            stds.append(r['c2_treatment_std'])
            pvals.append(r['c2_p_value'])

    baseline = results[0]['c2_control_mean']

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.errorbar(treatment_hours, means, yerr=stds, marker='o', markersize=10,
               linewidth=2.5, capsize=6, color='#2ca02c', label='Treatment')
    ax.axhline(baseline, color='gray', linestyle='--', linewidth=2, label='Control baseline')

    for i, (h, p) in enumerate(zip(treatment_hours, pvals)):
        if p < bonferroni_alpha:
            ax.plot(h, means[i], marker='*', markersize=20, color='gold', zorder=10)

    ax.set_xlabel('Treatment Duration (hours)', fontsize=12, fontweight='bold')
    ax.set_ylabel('C2 Scatter Index (a.u.)', fontsize=12, fontweight='bold')
    ax.set_title('Protein X Scatter Over Treatment Time', fontsize=14, fontweight='bold')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.set_xticks(treatment_hours)

    plt.tight_layout()
    plt.savefig(SVG_DIR / 'plot1_temporal_progression.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ plot1_temporal_progression.svg")
    plt.close()


def plot_effect_sizes():
    """Cohen's d effect sizes"""
    results, _ = load_results()

    treatment_labels = ['4h', '6h', '9h', '12h', '17h', '20h']
    effect_sizes = []
    for tl in treatment_labels:
        r = next((r for r in results if r['time_point'] == tl), None)
        if r:
            effect_sizes.append(r['c2_cohens_d'])

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ['#2ca02c' if d > 0 else '#d62728' for d in effect_sizes]
    bars = ax.bar(treatment_labels, effect_sizes, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)

    ax.axhline(0.2, color='gray', linestyle=':', label='Small (0.2)')
    ax.axhline(0.5, color='gray', linestyle='--', label='Medium (0.5)')
    ax.axhline(0.8, color='gray', linestyle='-', label='Large (0.8)')
    ax.axhline(0, color='black', linewidth=1)

    ax.set_xlabel('Treatment Duration', fontsize=12, fontweight='bold')
    ax.set_ylabel("Cohen's d", fontsize=12, fontweight='bold')
    ax.set_title('Effect Size by Treatment Duration', fontsize=14, fontweight='bold')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(SVG_DIR / 'plot2_effect_sizes.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ plot2_effect_sizes.svg")
    plt.close()


def plot_percent_change_waterfall():
    """Waterfall chart of percent changes"""
    results, bonferroni_alpha = load_results()

    treatment_labels = ['4h', '6h', '9h', '12h', '17h', '20h']
    changes, sigs = [], []
    for tl in treatment_labels:
        r = next((r for r in results if r['time_point'] == tl), None)
        if r:
            changes.append(r['c2_percent_change'])
            sigs.append(r['c2_p_value'] < bonferroni_alpha)

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ['#2ca02c' if s else '#98df8a' for s in sigs]
    bars = ax.bar(treatment_labels, changes, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)

    for i, (bar, sig) in enumerate(zip(bars, sigs)):
        if sig:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2, height + 0.5, '***',
                   ha='center', fontsize=14, fontweight='bold')

    ax.axhline(0, color='black', linewidth=1.5)
    ax.set_xlabel('Treatment Duration', fontsize=12, fontweight='bold')
    ax.set_ylabel('Percent Change from Controls (%)', fontsize=12, fontweight='bold')
    ax.set_title('Scatter Change from Baseline', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(SVG_DIR / 'plot6_percent_change_waterfall.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ plot6_percent_change_waterfall.svg")
    plt.close()


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

    exp1_data = load_from_experiment("experiments/outputs")
    exp2_data = load_from_experiment("experiments/outputs_exp2")
    exp3_data = load_from_experiment("experiments/outputs_exp3")
    return exp1_data, exp2_data, exp3_data


def map_condition_to_timepoint(condition_name):
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
    pooled = {}
    for exp_data in [exp1_data, exp2_data, exp3_data]:
        for condition_name, condition_data in exp_data.items():
            tp_name, tp_hours = map_condition_to_timepoint(condition_name)
            if tp_hours is None:
                continue
            if tp_name not in pooled:
                pooled[tp_name] = {'hours': tp_hours, 'c1': [], 'c2': []}
            pooled[tp_name]['c1'].extend(condition_data['c1'])
            pooled[tp_name]['c2'].extend(condition_data['c2'])
    return pooled


def plot_pvalue_heatmap():
    """P-value significance heatmap"""
    results, _ = load_results()

    fig, ax = plt.subplots(figsize=(10, 3))

    time_points = [r['time_point'] for r in results if not r['time_point'].startswith('recovery')]
    p_values = [r['c2_p_value'] for r in results if not r['time_point'].startswith('recovery')]

    sig_levels = []
    for p in p_values:
        if p < 0.001:
            sig_levels.append(4)
        elif p < 0.008333:
            sig_levels.append(3)
        elif p < 0.01:
            sig_levels.append(2)
        elif p < 0.05:
            sig_levels.append(1)
        else:
            sig_levels.append(0)

    data = np.array([sig_levels])
    im = ax.imshow(data, cmap='RdYlGn', aspect='auto', vmin=0, vmax=4)

    ax.set_xticks(range(len(time_points)))
    ax.set_xticklabels(time_points, fontsize=11)
    ax.set_yticks([0])
    ax.set_yticklabels(['Statistical\nSignificance'], fontsize=11)

    for i, (tp, p, sig) in enumerate(zip(time_points, p_values, sig_levels)):
        color = 'white' if sig >= 3 else 'black'
        if p < 0.001:
            text = 'p<0.001\n***'
        else:
            text = f'p={p:.4f}\n{"***" if sig >= 3 else ("**" if sig == 2 else ("*" if sig == 1 else "ns"))}'
        ax.text(i, 0, text, ha='center', va='center', fontsize=9, fontweight='bold', color=color)

    ax.set_title('Statistical Significance Across Treatment Durations\n(Bonferroni-corrected threshold: p < 0.008)',
                fontsize=13, fontweight='bold', pad=15)

    cbar = plt.colorbar(im, ax=ax, orientation='horizontal', pad=0.15, aspect=30)
    cbar.set_ticks([0, 1, 2, 3, 4])
    cbar.set_ticklabels(['ns\n(p≥0.05)', 'p<0.05', 'p<0.01', 'p<0.008\n(Bonf.)', 'p<0.001'])
    cbar.ax.tick_params(labelsize=8)

    plt.tight_layout()
    plt.savefig(SVG_DIR / 'plot3_pvalue_heatmap.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ plot3_pvalue_heatmap.svg")
    plt.close()


def plot_scatter_distributions(pooled_data):
    """Violin plots comparing distributions"""
    from matplotlib.patches import Patch

    fig, ax = plt.subplots(figsize=(14, 6))

    time_order = ['Controls', '4h', '6h', '9h', '12h', '17h', '20h']
    plot_data = []
    plot_labels = []
    plot_colors = []

    for tp in time_order:
        if tp in pooled_data:
            plot_data.append(pooled_data[tp]['c2'])
            plot_labels.append(f"{tp}\n(n={len(pooled_data[tp]['c2'])})")
            plot_colors.append('#3498DB' if tp == 'Controls' else '#E74C3C')

    parts = ax.violinplot(plot_data, positions=range(len(plot_data)), showmeans=True, showmedians=True)

    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(plot_colors[i])
        pc.set_alpha(0.7)
        pc.set_edgecolor('black')
        pc.set_linewidth(1)

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

    legend_elements = [
        Patch(facecolor='#3498DB', alpha=0.7, edgecolor='black', label='Controls'),
        Patch(facecolor='#E74C3C', alpha=0.7, edgecolor='black', label='Treatment')
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True, fancybox=True, shadow=True)

    plt.tight_layout()
    plt.savefig(SVG_DIR / 'plot4_scatter_distributions.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ plot4_scatter_distributions.svg")
    plt.close()


def plot_sample_sizes(pooled_data):
    """Sample sizes per condition"""
    results, _ = load_results()

    fig, ax = plt.subplots(figsize=(10, 5))

    time_points = ['Controls'] + [r['time_point'] for r in results if not r['time_point'].startswith('recovery')]
    sample_sizes = [len(pooled_data['Controls']['c2'])] + [r['n_treatment'] for r in results if not r['time_point'].startswith('recovery')]

    colors = ['#3498DB'] + ['#E74C3C'] * (len(time_points) - 1)

    bars = ax.bar(range(len(time_points)), sample_sizes, color=colors,
                  edgecolor='black', linewidth=1.5, alpha=0.8)

    for i, (bar, size) in enumerate(zip(bars, sample_sizes)):
        ax.text(i, size + 10, str(size), ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_xticks(range(len(time_points)))
    ax.set_xticklabels(time_points, fontsize=11)
    ax.set_ylabel('Number of Cells Analyzed', fontsize=12, fontweight='bold')
    ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
    ax.set_title('Statistical Power: Sample Sizes Across Conditions\n(Pooled from 3 Independent Experiments)',
                fontsize=13, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5, axis='y')

    total_cells = sum(sample_sizes)
    ax.text(0.98, 0.98, f'Total cells analyzed: {total_cells:,}',
           transform=ax.transAxes, ha='right', va='top', fontsize=11, fontweight='bold',
           bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    plt.tight_layout()
    plt.savefig(SVG_DIR / 'plot5_sample_sizes.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ plot5_sample_sizes.svg")
    plt.close()


SVG_CONTROLS_DIR = SVG_DIR / "plots_with_controls"


def plot_alt1_temporal_with_controls(pooled_data, results):
    """Alternative Plot 1: Temporal progression WITH controls at time=0"""
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(12, 6))

    all_conditions = ['Controls'] + ['4h', '6h', '9h', '12h', '17h', '20h']
    times, means, sems, colors, markers, significant = [], [], [], [], [], []

    for i, cond in enumerate(all_conditions):
        if cond in pooled_data:
            data = pooled_data[cond]
            times.append(data['hours'])
            means.append(np.mean(data['c2']))
            sems.append(np.std(data['c2']) / np.sqrt(len(data['c2'])))

            if cond == 'Controls':
                colors.append('#3498DB')
                markers.append('s')
                significant.append(False)
            else:
                result = next((r for r in results if r['time_point'] == cond), None)
                is_sig = result and result.get('c2_significant_bonferroni', False)
                significant.append(is_sig)
                colors.append('#E74C3C' if is_sig else '#95A5A6')
                markers.append('o')

    ax.plot(times, means, 'o-', linewidth=2.5, color='#34495E', markersize=0, alpha=0.3, zorder=1)

    for i, (time, mean, sem, color, marker, sig) in enumerate(zip(times, means, sems, colors, markers, significant)):
        ax.errorbar(time, mean, yerr=sem, fmt=marker, markersize=10, color=color,
                   markeredgecolor='black', markeredgewidth=2, capsize=5, capthick=2, elinewidth=2, zorder=2)
        if sig:
            ax.text(time, mean + sem + 0.5, '***', ha='center', va='bottom',
                   fontsize=14, fontweight='bold', color='#C0392B')
        if i == 0:
            ax.text(time, mean - sem - 0.8, 'Baseline\n(Controls)', ha='center', va='top',
                   fontsize=9, style='italic', color='#2C3E50',
                   bbox=dict(boxstyle='round,pad=0.4', facecolor='#ECF0F1', alpha=0.8))

    ax.set_xlabel('Time (hours)', fontsize=12, fontweight='bold')
    ax.set_ylabel('C2 Scatter Index\n(Protein X Dispersion)', fontsize=12, fontweight='bold')
    ax.set_title('Drug Effect on Protein Scattering: Time Course with Baseline\n(Pooled Analysis: 3 Experiments, 4,449 Cells)',
                fontsize=13, fontweight='bold', pad=15)

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
    ax.text(0.98, 0.02, '*** p < 0.008 (Bonferroni-corrected)',
           transform=ax.transAxes, ha='right', va='bottom',
           fontsize=9, style='italic', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(SVG_CONTROLS_DIR / 'plot_alt1_temporal_with_controls.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ plot_alt1_temporal_with_controls.svg")
    plt.close()


def plot_alt2_bar_chart_all_conditions(pooled_data, results):
    """Alternative Plot 2: Bar chart including controls"""
    fig, ax = plt.subplots(figsize=(12, 6))

    all_conditions = ['Controls', '4h', '6h', '9h', '12h', '17h', '20h']
    means, sems, colors, labels = [], [], [], []

    for cond in all_conditions:
        if cond in pooled_data:
            data = pooled_data[cond]
            means.append(np.mean(data['c2']))
            sems.append(np.std(data['c2']) / np.sqrt(len(data['c2'])))

            if cond == 'Controls':
                colors.append('#3498DB')
                labels.append(f'Controls\n(n={len(data["c2"])})')
            else:
                result = next((r for r in results if r['time_point'] == cond), None)
                is_sig = result and result.get('c2_significant_bonferroni', False)
                colors.append('#E74C3C' if is_sig else '#95A5A6')
                labels.append(f'{cond}\n(n={len(data["c2"])})')

    x_pos = np.arange(len(labels))
    bars = ax.bar(x_pos, means, yerr=sems, color=colors, edgecolor='black', linewidth=1.5,
                  alpha=0.8, capsize=5, error_kw={'elinewidth': 2, 'capthick': 2})

    for i, (bar, mean, sem) in enumerate(zip(bars, means, sems)):
        ax.text(i, mean + sem + 0.2, f'{mean:.2f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
        if i > 0:
            result = next((r for r in results if r['time_point'] == all_conditions[i]), None)
            if result and result.get('c2_significant_bonferroni', False):
                ax.text(i, mean + sem + 0.6, '***', ha='center', va='bottom',
                       fontsize=14, fontweight='bold', color='#C0392B')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel('C2 Scatter Index\n(Protein X Dispersion)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
    ax.set_title('Protein Scattering Across All Conditions\n(Mean ± SEM, *** = p<0.008)',
                fontsize=13, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5, axis='y')

    baseline = means[0]
    ax.axhline(baseline, color='#3498DB', linestyle='--', linewidth=2, alpha=0.5,
              label=f'Baseline (Controls: {baseline:.2f})')
    ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True)

    plt.tight_layout()
    plt.savefig(SVG_CONTROLS_DIR / 'plot_alt2_bar_chart_all_conditions.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ plot_alt2_bar_chart_all_conditions.svg")
    plt.close()


def plot_alt3_side_by_side_comparison(pooled_data, results):
    """Alternative Plot 3: Side-by-side comparison showing absolute values"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    all_conditions = ['Controls', '4h', '6h', '9h', '12h', '17h', '20h']
    means, sems, colors = [], [], []

    for cond in all_conditions:
        if cond in pooled_data:
            data = pooled_data[cond]
            means.append(np.mean(data['c2']))
            sems.append(np.std(data['c2']) / np.sqrt(len(data['c2'])))

            if cond == 'Controls':
                colors.append('#3498DB')
            else:
                result = next((r for r in results if r['time_point'] == cond), None)
                is_sig = result and result.get('c2_significant_bonferroni', False)
                colors.append('#E74C3C' if is_sig else '#95A5A6')

    x_pos = np.arange(len(all_conditions))

    # LEFT PLOT: Absolute values
    bars1 = ax1.bar(x_pos, means, yerr=sems, color=colors, edgecolor='black', linewidth=1.5,
                    alpha=0.8, capsize=5, error_kw={'elinewidth': 2})

    for i, (bar, mean, sem) in enumerate(zip(bars1, means, sems)):
        if i > 0:
            result = next((r for r in results if r['time_point'] == all_conditions[i]), None)
            if result and result.get('c2_significant_bonferroni', False):
                ax1.text(i, mean + sem + 0.2, '***', ha='center', va='bottom',
                        fontsize=12, fontweight='bold', color='#C0392B')

    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(all_conditions, fontsize=10)
    ax1.set_ylabel('C2 Scatter Index', fontsize=11, fontweight='bold')
    ax1.set_xlabel('Condition', fontsize=11, fontweight='bold')
    ax1.set_title('A) Absolute Scatter Values\n(Mean ± SEM)', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle=':', linewidth=0.5, axis='y')
    ax1.axhline(means[0], color='#3498DB', linestyle='--', linewidth=1.5, alpha=0.5)

    # RIGHT PLOT: Relative change from controls
    percent_changes = [0] + [((m - means[0]) / means[0]) * 100 for m in means[1:]]
    colors_right = ['#3498DB'] + colors[1:]
    bars2 = ax2.bar(x_pos, percent_changes, color=colors_right, edgecolor='black', linewidth=1.5, alpha=0.8)
    ax2.axhline(0, color='black', linestyle='-', linewidth=2)

    for i, (bar, pct) in enumerate(zip(bars2, percent_changes)):
        if i > 0:
            result = next((r for r in results if r['time_point'] == all_conditions[i]), None)
            if result and result.get('c2_significant_bonferroni', False):
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

    plt.suptitle('Protein Scattering: Absolute vs Relative Comparison\n(*** p < 0.008, Bonferroni-corrected)',
                fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout()
    plt.savefig(SVG_CONTROLS_DIR / 'plot_alt3_side_by_side_comparison.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ plot_alt3_side_by_side_comparison.svg")
    plt.close()


def plot_alt4_grouped_bar_with_stats(pooled_data, results):
    """Alternative Plot 4: Grouped bars with statistical annotations"""
    fig, ax = plt.subplots(figsize=(14, 7))

    all_conditions = ['Controls', '4h', '6h', '9h', '12h', '17h', '20h']
    x = np.arange(len(all_conditions))
    width = 0.6

    means, sems, p_values_text, colors = [], [], [], []

    for i, cond in enumerate(all_conditions):
        if cond in pooled_data:
            data = pooled_data[cond]
            means.append(np.mean(data['c2']))
            sems.append(np.std(data['c2']) / np.sqrt(len(data['c2'])))

            if cond == 'Controls':
                p_values_text.append('Baseline')
                colors.append('#3498DB')
            else:
                result = next((r for r in results if r['time_point'] == cond), None)
                is_sig = result and result.get('c2_significant_bonferroni', False)
                colors.append('#E74C3C' if is_sig else '#95A5A6')
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

    bars = ax.bar(x, means, width, yerr=sems, color=colors, edgecolor='black', linewidth=2,
                  alpha=0.85, capsize=6, error_kw={'elinewidth': 2.5, 'capthick': 2})

    for i, (bar, mean, sem, p_text) in enumerate(zip(bars, means, sems, p_values_text)):
        ax.text(i, mean + sem + 0.3, f'{mean:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
        ax.text(i, 9.5, p_text, ha='center', va='top', fontsize=8, style='italic',
               color='#2C3E50' if i == 0 else ('#C0392B' if '***' in p_text else '#7F8C8D'))
        n = len(pooled_data[all_conditions[i]]['c2'])
        ax.text(i, 9.2, f'n={n}', ha='center', va='top', fontsize=7, color='#34495E')

    ax.set_xticks(x)
    ax.set_xticklabels(all_conditions, fontsize=11, fontweight='bold')
    ax.set_ylabel('C2 Scatter Index (Protein X)', fontsize=12, fontweight='bold')
    ax.set_xlabel('Condition', fontsize=12, fontweight='bold')
    ax.set_title('Comprehensive Comparison: All Conditions with Statistics\n(Bars = Mean ± SEM, *** = Bonferroni-corrected p<0.008)',
                fontsize=13, fontweight='bold', pad=20)

    baseline = means[0]
    ax.axhline(baseline, color='#3498DB', linestyle='--', linewidth=2, alpha=0.4,
              label=f'Controls Baseline: {baseline:.2f}')
    ax.grid(True, alpha=0.2, linestyle=':', linewidth=0.5, axis='y')
    ax.set_ylim(9, 15)
    ax.legend(loc='upper left', fontsize=10, frameon=True, fancybox=True, shadow=True)
    ax.axhspan(baseline - sems[0], baseline + sems[0], color='#3498DB', alpha=0.1, zorder=0)

    plt.tight_layout()
    plt.savefig(SVG_CONTROLS_DIR / 'plot_alt4_grouped_bar_with_stats.svg', format='svg', bbox_inches='tight')
    print(f"  ✅ plot_alt4_grouped_bar_with_stats.svg")
    plt.close()


def main():
    """Generate all plots in SVG format"""
    print("\n" + "="*70)
    print("CONVERTING ALL PLOTS TO SVG (VECTOR FORMAT)")
    print("="*70)
    print("\nSVG files can be edited in:")
    print("  - Adobe Illustrator")
    print("  - Inkscape (free)")
    print("  - Figma")
    print("  - Microsoft PowerPoint (insert as image, then ungroup)")
    print("  - Affinity Designer")
    print()

    # Create directories
    SVG_DIR.mkdir(parents=True, exist_ok=True)
    SVG_DUAL_DIR.mkdir(parents=True, exist_ok=True)
    SVG_CONTROLS_DIR.mkdir(parents=True, exist_ok=True)

    print("Generating dual-channel plots (SVG)...")
    plot_dual_channel_temporal()
    plot_percent_change_comparison()
    plot_scatter_plot_c1_vs_c2()
    plot_effect_size_comparison()
    plot_reversibility_analysis()

    print("\nGenerating single-channel plots (SVG)...")
    plot_temporal_progression()
    plot_effect_sizes()
    plot_pvalue_heatmap()
    plot_percent_change_waterfall()

    print("\nGenerating distribution/sample plots (SVG)...")
    try:
        exp1_data, exp2_data, exp3_data = load_pooled_data()
        pooled_data = pool_by_timepoint(exp1_data, exp2_data, exp3_data)
        plot_scatter_distributions(pooled_data)
        plot_sample_sizes(pooled_data)
    except Exception as e:
        print(f"  ⚠️ Could not generate distribution plots: {e}")

    print("\nGenerating plots with controls (SVG)...")
    try:
        with open('results/pooled_analysis_results.json', 'r') as f:
            analysis_results = json.load(f)
        results = analysis_results['results']

        exp1_data, exp2_data, exp3_data = load_pooled_data()
        pooled_data = pool_by_timepoint(exp1_data, exp2_data, exp3_data)

        plot_alt1_temporal_with_controls(pooled_data, results)
        plot_alt2_bar_chart_all_conditions(pooled_data, results)
        plot_alt3_side_by_side_comparison(pooled_data, results)
        plot_alt4_grouped_bar_with_stats(pooled_data, results)
    except Exception as e:
        print(f"  ⚠️ Could not generate plots with controls: {e}")

    print("\n" + "="*70)
    print(f"✅ ALL SVG FILES SAVED TO: {SVG_DIR.resolve()}")
    print("="*70)

    # List all generated files
    svg_files = list(SVG_DIR.rglob("*.svg"))
    print(f"\nGenerated {len(svg_files)} SVG files:")
    for f in sorted(svg_files):
        print(f"  • {f.relative_to(SVG_DIR.parent.parent)}")


if __name__ == '__main__':
    main()
