#!/usr/bin/env python3
"""
Temporal progression, dual channel, gene-named styling.

Treatment course shown as 4h -> 6h -> 9h -> 20h (12h and 17h excluded), running
straight on into the single 20h-treatment recovery point after washout (the 12h
and 17h recovery conditions are excluded too).

Significance convention matches the earlier published version of this figure:
*** is drawn where p < 0.01 AND the mean is above the channel's control baseline
(i.e. significant scattering increases only).
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt

RESULTS_JSON = Path('results/pooled_analysis_results.json')
OUT_PNG = Path('results/plots/png/dual_channel')
OUT_SVG = Path('results/plots/svg/dual_channel')

# c2 = HTT-mNG (green), c1 = GOLGA2 (red)
GREEN = '#4a9d4a'
RED = '#cf5b52'

TREATMENT = [('4h', '4h'), ('6h', '6h'), ('9h', '9h'), ('20h', '20h')]
# Only the 20h-treatment recovery condition is shown, so the curve reads as one
# continuous course: 20h of drug, then washout and recovery.
RECOVERY = [('recovery_20h', '20h recovery')]

SIG_P = 0.01


def load_results():
    with open(RESULTS_JSON) as f:
        data = json.load(f)
    return {r['time_point']: r for r in data['results']}, data['results'][0]


def series(rows, keys, channel):
    means = [rows[k][f'{channel}_treatment_mean'] for k, _ in keys]
    pvals = [rows[k][f'{channel}_p_value'] for k, _ in keys]
    return means, pvals


def main():
    rows, first = load_results()
    OUT_PNG.mkdir(parents=True, exist_ok=True)
    OUT_SVG.mkdir(parents=True, exist_ok=True)

    c1_base = first['c1_control_mean']
    c2_base = first['c2_control_mean']

    # One continuous axis: recovery follows straight on from the 20h point
    t_x = list(range(len(TREATMENT)))
    r_x = list(range(len(TREATMENT), len(TREATMENT) + len(RECOVERY)))

    c1_t, c1_tp = series(rows, TREATMENT, 'c1')
    c2_t, c2_tp = series(rows, TREATMENT, 'c2')
    c1_r, c1_rp = series(rows, RECOVERY, 'c1')
    c2_r, c2_rp = series(rows, RECOVERY, 'c2')

    fig, ax = plt.subplots(figsize=(11, 5.5))

    ax.axhline(c2_base, color=GREEN, linestyle=':', linewidth=1,
               alpha=0.7, label='HTT-mNG baseline')
    ax.axhline(c1_base, color=RED, linestyle=':', linewidth=1,
               alpha=0.7, label='GOLGA2 baseline')

    ax.plot(t_x, c2_t, marker='s', markersize=7, linewidth=2,
            color=GREEN, label='HTT-mNG')
    ax.plot(t_x, c1_t, marker='o', markersize=7, linewidth=2,
            color=RED, label='GOLGA2')

    join_x = [t_x[-1]] + r_x
    ax.plot(join_x, [c2_t[-1]] + c2_r, marker='s', markersize=7, linewidth=2,
            color=GREEN, linestyle='--', alpha=0.85)
    ax.plot(join_x, [c1_t[-1]] + c1_r, marker='o', markersize=7, linewidth=2,
            color=RED, linestyle='--', alpha=0.85)

    def stars(xs, means, pvals, baseline, color, offset):
        for x, m, p in zip(xs, means, pvals):
            if p < SIG_P and m > baseline:
                ax.text(x, m + offset, '***', color=color, ha='center',
                        va='bottom', fontsize=11, fontweight='bold')

    stars(t_x, c2_t, c2_tp, c2_base, GREEN, 0.30)
    stars(t_x, c1_t, c1_tp, c1_base, RED, 0.85)
    stars(r_x, c2_r, c2_rp, c2_base, GREEN, 0.30)
    stars(r_x, c1_r, c1_rp, c1_base, RED, 0.85)

    # Shade the recovery phase rather than breaking the axis
    ax.axvspan(t_x[-1], r_x[-1] + 0.45, color='#bbbbbb', alpha=0.10, zorder=0)
    ax.text(sum(t_x) / len(t_x), 17.6, 'Drug treatment', ha='center', va='center',
            fontsize=10, color='#555555', style='italic')
    ax.text(sum(r_x) / len(r_x), 17.6, 'Recovery', ha='center',
            va='center', fontsize=10, color='#555555', style='italic')

    ax.set_xticks(t_x + r_x)
    ax.set_xticklabels([lbl for _, lbl in TREATMENT] + [lbl for _, lbl in RECOVERY])
    ax.set_xlim(-0.45, r_x[-1] + 0.45)
    ax.set_ylim(7.5, 18.5)
    ax.set_yticks([8, 10, 12, 14, 16, 18])
    ax.set_xlabel('Treatment duration / recovery after washout', fontsize=11)
    ax.set_ylabel('Scatter Index (a.u.)', fontsize=11)

    ax.grid(True, axis='y', alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color('#333333')
        spine.set_linewidth(0.8)

    ax.legend(loc='upper left', frameon=False, fontsize=10)

    plt.tight_layout()
    plt.savefig(OUT_PNG / 'temporal_progression_dual.png', dpi=300, bbox_inches='tight')
    plt.savefig(OUT_SVG / 'temporal_progression_dual.svg', format='svg', bbox_inches='tight')
    print(f"Saved: {OUT_PNG / 'temporal_progression_dual.png'}")
    print(f"Saved: {OUT_SVG / 'temporal_progression_dual.svg'}")
    plt.close()


if __name__ == '__main__':
    main()
