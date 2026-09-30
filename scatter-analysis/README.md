# Dual-Channel Protein Scatter Analysis

Quantitative analysis of drug-induced protein scattering in fluorescent microscopy images, with multi-experiment pooling and recovery analysis.

**Associated paper**: *[Paper title / DOI to be added]*

## Overview

This pipeline analyzes how 2-bromopalmitate (2-BP) treatment affects the spatial distribution of two fluorescent proteins in cells:

- **C1 (Red)** -- Structure Y (reference structure)
- **C2 (Green)** -- Protein X (primary target)
- **C3 (Blue)** -- Positioning marker (segmentation only)

The analysis measures *scatter index* (intensity-weighted standard deviation from cell center) across 10 temporal conditions: 1 control, 6 treatment time points (4h--20h), and 3 recovery conditions (12h, 17h, 20h washout).

### Key Findings

| Finding | Verdict | Evidence |
|---------|---------|----------|
| Drug increases Protein X scattering | Strongly supported | 4/6 time points significant (peak: +15.0% at 9h, d=0.44) |
| C1 and C2 co-scatter at peak times | Supported | Both channels significant at 9h and 20h |
| Drug effect is reversible | Strongly supported | Both channels recover to baseline by 20h washout |

- **Sample size**: 4,449 cells pooled across 3 independent experiments
- **Statistical method**: Welch's t-test with Bonferroni correction (alpha = 0.005556 for 9 comparisons)
- **Effect sizes**: Cohen's d reported for all comparisons


## Repository Structure

```
orchid/
├── pipeline/                              # 7-stage image processing pipeline
│   ├── stages/
│   │   ├── stage_01_validate.py               # Data loading and validation
│   │   ├── stage_02_segment.py                # Cell segmentation (contour detection)
│   │   ├── stage_03_register.py               # Cross-channel cell registration
│   │   ├── stage_04_detect_particles.py       # Particle detection (Otsu)
│   │   ├── stage_05_calculate_metrics.py      # Scatter index and correlation metrics
│   │   ├── stage_06_analyze.py                # Statistical comparison
│   │   └── stage_07_visualize.py              # Report generation
│   └── cli/                               # CLI commands for run_pipeline.py
├── run_pipeline.py                        # CLI entry point (single experiment)
├── process_experiment.py                  # Run stages 1-5 on one experiment
├── analyze_pooled_experiments.py          # Pool 3 experiments, run statistical tests
├── create_pooled_plots.py                 # Single-channel temporal plots (6 figures)
├── create_dual_channel_plots.py           # Dual-channel and recovery plots (5 figures)
├── create_plots_with_controls.py          # Alternative plots with control data points
├── convert_plots_to_pdf.py               # PNG to PDF conversion
├── convert_plots_to_svg.py               # Regenerate all plots as editable SVG
├── context/
│   └── Experiment Analysis Summary.txt    # Experiment protocol description
├── requirements.txt
└── README.md
```

**Note**: The `experiments/` directory containing raw microscopy data (~1.2 GB) is not included in this repository. See [Data Availability](#data-availability).

## Reproducing the Analysis

### Prerequisites

- Python 3.10+
- Raw experiment data placed in `experiments/` (see [Data Availability](#data-availability))

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Process each experiment (stages 1-5)

From the repository root:

```bash
python process_experiment.py "experiments/Experiment 1" experiments/outputs
python process_experiment.py "experiments/Experiment 2" experiments/outputs_exp2
python process_experiment.py "experiments/Experiment 3" experiments/outputs_exp3
```

Each run processes ~20 samples per condition through cell segmentation, channel registration, particle detection, and scatter metric calculation. Boundary cells (within 10% of image edge) are excluded.

### 3. Run pooled statistical analysis

```bash
python analyze_pooled_experiments.py
```

Pools scatter indices across all 3 experiments by time point, runs Welch's t-tests (treatment vs. control), applies Bonferroni correction, and writes results to `results/pooled_analysis_results.json`. Creates the `results/` directory if needed.

### 4. Generate figures

```bash
python create_pooled_plots.py          # 6 single-channel C2 plots
python create_dual_channel_plots.py    # 5 dual-channel + recovery plots
python create_plots_with_controls.py   # Alternative plots with control data points
```

### 5. Convert to publication formats (optional)

```bash
python convert_plots_to_pdf.py         # PNG -> PDF
python convert_plots_to_svg.py         # Editable SVG (text preserved as text)
```

## Methods

### Image Processing Pipeline (Stages 1-5)

1. **Validation** -- Verifies 3 channels (C1, C2, C3) present per sample
2. **Segmentation** -- Detects cells using contour detection on C3; excludes boundary cells (10% edge margin)
3. **Registration** -- Maps segmented cells across channels using shared bounding boxes
4. **Particle Detection** -- Identifies fluorescent particles using Otsu thresholding
5. **Metrics** -- Computes scatter index (intensity-weighted standard deviation from cell centroid), pixel Pearson correlation, and Manders M1/M2 colocalization coefficients

### Statistical Analysis

- **Comparison**: Each of 9 conditions (6 treatment + 3 recovery) vs. pooled controls
- **Test**: Welch's t-test (does not assume equal variance)
- **Multiple comparison correction**: Bonferroni (alpha = 0.05 / 9 = 0.005556)
- **Effect size**: Cohen's d with pooled standard deviation

### Known Limitations

- **Pseudo-replication**: Cells from the same sample are treated as independent observations; a mixed-effects model would better account for the nested structure
- **No cell-level correlation**: Co-scattering is inferred from population-level channel statistics, not from paired C1/C2 measurements within individual cells
- **Cross-sectional design**: Different cells are measured at each time point (not longitudinal tracking)

## Data Availability

Raw microscopy data is available from *[contact / repository to be added]*. Place the experiment directories as:

```
experiments/
├── Experiment 1/    # Controls + 6 treatment + 3 recovery conditions
├── Experiment 2/    # Independent replicate
└── Experiment 3/    # Independent replicate
```

Each experiment contains ~20 TIFF samples per condition, with 3 channels (C1, C2, C3) per sample.

## License

*[License to be added]*
