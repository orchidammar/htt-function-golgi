#!/usr/bin/env python3
"""
Pool scatter index data across all 3 experiments and analyze Treatment vs Controls
for each time point.

Strategy:
- For each time point (4h, 6h, 9h, 12h, 17h, 20h):
  * Pool Controls from Exp1 + Exp2 + Exp3
  * Pool Treatment from Exp1 + Exp2 + Exp3
  * Compare: Treatment vs Controls (t-test)
  * Calculate: Cohen's d effect size
- Apply Bonferroni correction for multiple comparisons
"""

import json
import numpy as np
from pathlib import Path
from scipy import stats as scipy_stats

def load_scatter_indices(output_dir):
    """
    Load all scatter index data from an experiment's output directory.
    Returns dict: {condition: {sample: {cell_id: {c1_scatter, c2_scatter}}}}
    """
    metrics_dir = Path(output_dir) / "05_metrics"
    if not metrics_dir.exists():
        return {}

    data = {}

    # Iterate through condition directories
    for condition_dir in metrics_dir.iterdir():
        if not condition_dir.is_dir():
            continue

        condition_name = condition_dir.name
        data[condition_name] = {}

        # Iterate through sample directories
        for sample_dir in condition_dir.iterdir():
            if not sample_dir.is_dir():
                continue

            sample_name = sample_dir.name
            correlation_dir = sample_dir / "correlation"

            if not correlation_dir.exists():
                continue

            data[condition_name][sample_name] = {}

            # Load correlation JSON files (one per cell)
            for corr_file in correlation_dir.glob("*_correlation.json"):
                with open(corr_file, 'r') as f:
                    corr_data = json.load(f)

                cell_id = corr_file.stem.replace('_correlation', '')

                # Extract scatter indices
                data[condition_name][sample_name][cell_id] = {
                    'c1_scatter': corr_data.get('c1_scatter_index'),
                    'c2_scatter': corr_data.get('c2_scatter_index')
                }

    return data

def pool_by_time_point(exp1_data, exp2_data, exp3_data):
    """
    Pool scatter data by time point across all experiments.
    Returns dict: {time_point: {'controls': [values], 'treatment': [values]}}
    """
    # Time point mapping (normalize condition names)
    time_point_map = {
        'Controls': 'controls',
        '4 hours treatment': '4h',
        '6 hours treatment': '6h',
        '9 hours treatment': '9h',
        '12 hours treatment': '12h',
        '17 hours treatment': '17h',
        '20 hours treatment': '20h',
        'Recovery after 12h treatment': 'recovery_12h',
        'Recovery after 12 hours': 'recovery_12h',
        'Recovery after 17h treatment': 'recovery_17h',
        'Recovery after 17 hours': 'recovery_17h',
        'Recovery after 20h treatment': 'recovery_20h',
        'Recovery after 20 hours': 'recovery_20h'
    }

    pooled = {}

    # Process each experiment
    for exp_data, exp_name in [(exp1_data, 'Exp1'), (exp2_data, 'Exp2'), (exp3_data, 'Exp3')]:
        for condition_name, condition_data in exp_data.items():
            # Map to standardized time point
            time_point = time_point_map.get(condition_name)
            if not time_point:
                print(f"Warning: Unknown condition '{condition_name}' in {exp_name}")
                continue

            # Initialize time point if needed
            if time_point not in pooled:
                pooled[time_point] = {'c1_scatter': [], 'c2_scatter': []}

            # Collect all scatter values from all cells in all samples
            for sample_name, sample_data in condition_data.items():
                for cell_id, cell_data in sample_data.items():
                    if cell_data['c1_scatter'] is not None:
                        pooled[time_point]['c1_scatter'].append(cell_data['c1_scatter'])
                    if cell_data['c2_scatter'] is not None:
                        pooled[time_point]['c2_scatter'].append(cell_data['c2_scatter'])

    return pooled

def compare_treatment_vs_controls(pooled_data):
    """
    For each treatment time point, compare against Controls.
    Returns comparison results with p-values and effect sizes.
    """
    controls_c1 = np.array(pooled_data.get('controls', {}).get('c1_scatter', []))
    controls_c2 = np.array(pooled_data.get('controls', {}).get('c2_scatter', []))

    print(f"\nBaseline (Controls):")
    print(f"  C1 scatter: {np.mean(controls_c1):.2f} ± {np.std(controls_c1):.2f} (n={len(controls_c1)})")
    print(f"  C2 scatter: {np.mean(controls_c2):.2f} ± {np.std(controls_c2):.2f} (n={len(controls_c2)})")

    # Treatment time points AND recovery conditions to compare
    treatment_times = ['4h', '6h', '9h', '12h', '17h', '20h',
                       'recovery_12h', 'recovery_17h', 'recovery_20h']

    results = []

    for time_point in treatment_times:
        if time_point not in pooled_data:
            print(f"\nWarning: No data for {time_point}")
            continue

        treatment_c1 = np.array(pooled_data[time_point]['c1_scatter'])
        treatment_c2 = np.array(pooled_data[time_point]['c2_scatter'])

        # Prepare result dict
        result = {
            'time_point': time_point,
            'n_control': len(controls_c2),
            'n_treatment': len(treatment_c2)
        }

        # C1 scatter comparison (Structure Y - RED channel)
        if len(treatment_c1) > 0 and len(controls_c1) > 0:
            t_stat_c1, p_value_c1 = scipy_stats.ttest_ind(treatment_c1, controls_c1)

            # Cohen's d for C1
            pooled_std_c1 = np.sqrt(((len(controls_c1) - 1) * np.std(controls_c1)**2 +
                                      (len(treatment_c1) - 1) * np.std(treatment_c1)**2) /
                                     (len(controls_c1) + len(treatment_c1) - 2))
            cohens_d_c1 = (np.mean(treatment_c1) - np.mean(controls_c1)) / pooled_std_c1 if pooled_std_c1 > 0 else 0

            percent_change_c1 = ((np.mean(treatment_c1) - np.mean(controls_c1)) / np.mean(controls_c1)) * 100

            result.update({
                'c1_control_mean': float(np.mean(controls_c1)),
                'c1_control_std': float(np.std(controls_c1)),
                'c1_treatment_mean': float(np.mean(treatment_c1)),
                'c1_treatment_std': float(np.std(treatment_c1)),
                'c1_percent_change': percent_change_c1,
                'c1_p_value': float(p_value_c1),
                'c1_cohens_d': float(cohens_d_c1)
            })

        # C2 scatter comparison (Protein X - GREEN channel)
        if len(treatment_c2) > 0 and len(controls_c2) > 0:
            t_stat_c2, p_value_c2 = scipy_stats.ttest_ind(treatment_c2, controls_c2)

            # Cohen's d for C2
            pooled_std_c2 = np.sqrt(((len(controls_c2) - 1) * np.std(controls_c2)**2 +
                                      (len(treatment_c2) - 1) * np.std(treatment_c2)**2) /
                                     (len(controls_c2) + len(treatment_c2) - 2))
            cohens_d_c2 = (np.mean(treatment_c2) - np.mean(controls_c2)) / pooled_std_c2 if pooled_std_c2 > 0 else 0

            percent_change_c2 = ((np.mean(treatment_c2) - np.mean(controls_c2)) / np.mean(controls_c2)) * 100

            result.update({
                'c2_control_mean': float(np.mean(controls_c2)),
                'c2_control_std': float(np.std(controls_c2)),
                'c2_treatment_mean': float(np.mean(treatment_c2)),
                'c2_treatment_std': float(np.std(treatment_c2)),
                'c2_percent_change': percent_change_c2,
                'c2_p_value': float(p_value_c2),
                'c2_cohens_d': float(cohens_d_c2)
            })

        results.append(result)

    return results

def apply_bonferroni_correction(results):
    """Apply Bonferroni correction to p-values"""
    n_comparisons = len(results)
    alpha = 0.05
    bonferroni_alpha = alpha / n_comparisons

    print(f"\nBonferroni correction:")
    print(f"  Number of comparisons: {n_comparisons}")
    print(f"  Corrected alpha: {bonferroni_alpha:.6f}")

    for result in results:
        result['bonferroni_alpha'] = bonferroni_alpha
        if 'c1_p_value' in result:
            result['c1_significant_bonferroni'] = result['c1_p_value'] < bonferroni_alpha
        if 'c2_p_value' in result:
            result['c2_significant_bonferroni'] = result['c2_p_value'] < bonferroni_alpha

    return results, bonferroni_alpha

def print_results_table(results, bonferroni_alpha):
    """Print formatted results table"""
    print("\n" + "="*120)
    print("POOLED ANALYSIS: TREATMENT VS CONTROLS (All 3 Experiments Combined)")
    print("="*120)
    print(f"Statistical threshold: p < {bonferroni_alpha:.6f} (Bonferroni-corrected)")
    print()

    # C1 (Structure Y - RED) results
    print("\n=== C1 SCATTER (Structure Y - RED channel) ===")
    print(f"{'Time':6} {'Ctrl Mean':>10} {'Trt Mean':>10} {'Change':>10} {'p-value':>12} {'d':>7} {'Sig':>5} {'n':>6}")
    print("-"*70)
    for r in results:
        if 'c1_p_value' in r:
            sig = '***' if r.get('c1_significant_bonferroni', False) else ('**' if r['c1_p_value'] < 0.01 else ('*' if r['c1_p_value'] < 0.05 else 'ns'))
            print(f"{r['time_point']:6} "
                  f"{r['c1_control_mean']:10.2f} "
                  f"{r['c1_treatment_mean']:10.2f} "
                  f"{r['c1_percent_change']:+9.1f}% "
                  f"{r['c1_p_value']:12.6f} "
                  f"{r['c1_cohens_d']:7.3f} "
                  f"{sig:>5} "
                  f"{r['n_treatment']:>6}")

    # C2 (Protein X - GREEN) results
    print("\n=== C2 SCATTER (Protein X - GREEN channel) ===")
    print(f"{'Time':6} {'Ctrl Mean':>10} {'Trt Mean':>10} {'Change':>10} {'p-value':>12} {'d':>7} {'Sig':>5} {'n':>6}")
    print("-"*70)
    for r in results:
        if 'c2_p_value' in r:
            sig = '***' if r.get('c2_significant_bonferroni', False) else ('**' if r['c2_p_value'] < 0.01 else ('*' if r['c2_p_value'] < 0.05 else 'ns'))
            print(f"{r['time_point']:6} "
                  f"{r['c2_control_mean']:10.2f} "
                  f"{r['c2_treatment_mean']:10.2f} "
                  f"{r['c2_percent_change']:+9.1f}% "
                  f"{r['c2_p_value']:12.6f} "
                  f"{r['c2_cohens_d']:7.3f} "
                  f"{sig:>5} "
                  f"{r['n_treatment']:>6}")

    print()
    print("="*120)
    print("VERDICT:")
    print("="*120)

    # C1 verdict
    c1_significant = [r for r in results if r.get('c1_significant_bonferroni', False)]
    print("\n📍 STRUCTURE Y (C1 - RED channel):")
    if c1_significant:
        print("✅ Significant scatter increases:")
        for r in c1_significant:
            print(f"   • {r['time_point']:5}: {r['c1_percent_change']:+.1f}% "
                  f"(p={r['c1_p_value']:.6f}, d={r['c1_cohens_d']:.2f})")
    else:
        print("❌ No significant scatter increases at Bonferroni-corrected level")

    # C2 verdict
    c2_significant = [r for r in results if r.get('c2_significant_bonferroni', False)]
    print("\n📍 PROTEIN X (C2 - GREEN channel):")
    if c2_significant:
        print("✅ Significant scatter increases:")
        for r in c2_significant:
            print(f"   • {r['time_point']:5}: {r['c2_percent_change']:+.1f}% "
                  f"(p={r['c2_p_value']:.6f}, d={r['c2_cohens_d']:.2f})")
    else:
        print("❌ No significant scatter increases at Bonferroni-corrected level")

    # Combined verdict
    both_significant = [r for r in results if r.get('c1_significant_bonferroni', False) and r.get('c2_significant_bonferroni', False)]
    print("\n📍 CO-SCATTERING (Both channels together):")
    if both_significant:
        print("✅ Time points where BOTH channels scatter significantly:")
        for r in both_significant:
            print(f"   • {r['time_point']:5}: C1 {r['c1_percent_change']:+.1f}%, C2 {r['c2_percent_change']:+.1f}%")
    else:
        print("❌ No time points show significant scatter in BOTH channels simultaneously")

    print("="*120)

def main():
    print("\n" + "="*100)
    print("LOADING DATA FROM ALL 3 EXPERIMENTS")
    print("="*100)

    # Load data from each experiment
    print("\nLoading Experiment 1...")
    exp1_data = load_scatter_indices("experiments/outputs")
    print(f"  Found {len(exp1_data)} conditions")

    print("\nLoading Experiment 2...")
    exp2_data = load_scatter_indices("experiments/outputs_exp2")
    print(f"  Found {len(exp2_data)} conditions")

    print("\nLoading Experiment 3...")
    exp3_data = load_scatter_indices("experiments/outputs_exp3")
    print(f"  Found {len(exp3_data)} conditions")

    # Pool by time point
    print("\n" + "="*100)
    print("POOLING DATA BY TIME POINT")
    print("="*100)
    pooled_data = pool_by_time_point(exp1_data, exp2_data, exp3_data)

    for time_point, data in sorted(pooled_data.items()):
        print(f"  {time_point:15}: {len(data['c2_scatter']):4} cells (C2)")

    # Compare treatment vs controls
    print("\n" + "="*100)
    print("STATISTICAL COMPARISONS")
    print("="*100)
    results = compare_treatment_vs_controls(pooled_data)

    # Apply Bonferroni correction
    results, bonferroni_alpha = apply_bonferroni_correction(results)

    # Print results table
    print_results_table(results, bonferroni_alpha)

    # Save results
    output_file = Path("results/pooled_analysis_results.json")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump({
            'results': results,
            'bonferroni_alpha': bonferroni_alpha,
            'pooled_sample_sizes': {tp: len(data['c2_scatter']) for tp, data in pooled_data.items()}
        }, f, indent=2)

    print(f"\n✅ Results saved to: {output_file}")

if __name__ == '__main__':
    main()
