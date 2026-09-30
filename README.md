# htt-function-golgi

Post-processing and outlier screening of per-cell Golgi morphology
measurements from segmented microscopy images.

The script `analyse_sub_folders.py` walks through a directory tree of
experiment folders, collects the measurement CSV files in each folder, and
merges them into one row per cell. For each cell it computes three Golgi
morphology features:

- how large the Golgi convex hull is relative to the nucleus,
- how compactly the Golgi fills its convex hull,
- how many Golgi fragments (stacks) the cell contains.

Each feature is then standardised with a z-score, and cells that deviate
by more than 2 standard deviations on any feature are flagged as
out-of-range. The analysis runs once per experiment folder (batch) and then
again on all batches pooled, separately for each experimental channel.

Segmentation and measurement are performed beforehand. This script does not
process images. It only reads the resulting per-object measurements.

## Channels and batches

Each experiment folder is assigned to an experimental channel based on its
folder path, which must contain one of the following names
(case-insensitive):

- `KO`
- `HAP40KO`
- `Q23`
- `Q145`

Folders are numbered as batches independently for each channel, in the
order they are found. For example, the first `KO` folder becomes batch
`KO_0`, the second `KO_1`, and the first `Q23` folder `Q23_0`.

## Input data

Each experiment folder must contain the following CSV files:

| File                      | Required columns                                     | Used for                     |
|---------------------------|------------------------------------------------------|------------------------------|
| `cell_outline.csv`        | `ImageNumber`, `ObjectNumber`                        | Defines the cells            |
| `Golgi_convex_hull.csv`   | `ImageNumber`, `ObjectNumber`, `AreaShape_Area`      | Golgi convex hull area       |
| `nuclei.csv`              | `ImageNumber`, `ObjectNumber`, `AreaShape_Area`      | Nucleus area                 |
| `golgi_disconnected.csv`  | `ImageNumber`, `ObjectNumber`, `AreaShape_Area`      | Golgi stack area             |
| `golgi_stack.csv`         | `ImageNumber`, `ObjectNumber`, `Parent_cell_outline` | Golgi fragment count         |

Objects are matched to cells by `ImageNumber` and `ObjectNumber`. Starting
from `cell_outline.csv`, the other tables are left-joined, so every cell is
kept. In `golgi_stack.csv`, `Parent_cell_outline` gives the `ObjectNumber`
of the cell each Golgi fragment belongs to. The convex hull area is read
from the input measurements and is not computed by the script. Any other
CSV files in the folders are ignored.

An example layout:

```
data/
├── KO/
│   ├── experiment_1/
│   │   ├── cell_outline.csv
│   │   ├── Golgi_convex_hull.csv
│   │   ├── nuclei.csv
│   │   ├── golgi_disconnected.csv
│   │   └── golgi_stack.csv
│   └── experiment_2/
│       └── ...
├── HAP40KO/
│   └── ...
├── Q23/
│   └── ...
└── Q145/
    └── ...
```

## Features

For each cell, the following measurements are taken from the input data:

| Symbol                | Measurement                                 | Output column       |
|-----------------------|---------------------------------------------|---------------------|
| $A_{hull}$            | Golgi convex hull area                      | `convex_hull_area`  |
| $A_{nucleus}$         | Nucleus area                                | `nucleus_area`      |
| $A_{stack}$           | Golgi stack area                            | `golgi_stack_area`  |
| $N_{fragments}$       | Number of Golgi fragments in the cell       | `fragments_count`   |

From these, three features are derived.

**Golgi-to-nucleus area ratio** $R$ (`golgi_area`), which measures Golgi
size normalised to nucleus size:

$$R = \frac{A_{hull}}{A_{nucleus}}$$

**Golgi compactness** $C$ (`golgi_compactness`), the fraction of the convex
hull that is filled by Golgi. A value close to 1 indicates a compact Golgi,
and lower values indicate a dispersed or fragmented Golgi:

$$C = \frac{A_{stack}}{A_{hull}}$$

**Fragment count** $N_{fragments}$ (`fragments_count`), the number of Golgi
objects in `golgi_stack.csv` assigned to the cell. Higher values indicate
stronger Golgi fragmentation.

## Outlier screening

Each feature $x$ is converted to a z-score:

$$z = \frac{x - \mu}{\sigma}$$

where $\mu$ is the mean and $\sigma$ the sample standard deviation of that
feature across the current analysis population. The population is either
all cells of one experiment folder, or all pooled cells of one channel.
This produces the columns `golgi_area_zscore`, `golgi_compactness_zscore`,
and `fragments_count_zscore`.

A cell is classified as **out-of-range** (`in_range = "no"`) if the absolute
z-score of **any** of the three features exceeds 2:

$$|z_R| > 2 \quad \lor \quad |z_C| > 2 \quad \lor \quad |z_{N_{fragments}}| > 2$$

Otherwise it is **in-range** (`in_range = "yes"`).

Outliers are defined statistically, relative to the distribution of the
analysed population, not by fixed biological thresholds. "Out-of-range"
means a cell is unusual within its batch or channel. It does not
necessarily mean the cell is biologically abnormal.

## Installation

Requires Python 3.9 or newer.

```bash
git clone <repository-url>
cd htt-function-golgi
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

The script searches the **current working directory** and all of its
subdirectories, so run it from the top-level data folder:

```bash
cd /path/to/data
python /path/to/analyse_sub_folders.py
```

## Output

The analysis runs at two levels, and each level writes three CSV files.

**Per batch** (written into each experiment folder). Z-scores are relative
to the cells of that batch.

| File                        | Contents                                  |
|-----------------------------|-------------------------------------------|
| `analysis_with_zscore.csv`  | All cells with features and z-scores      |
| `analysis_in_range.csv`     | In-range cells only                       |
| `analysis_out_of_range.csv` | Out-of-range cells only                   |

**Per channel, pooled across batches** (written into the directory the
script is run from, for each of `KO`, `HAP40KO`, `Q23` and `Q145`). Z-scores
are relative to all cells of that channel, so channels are never mixed when
calculating the mean and standard deviation.

| File                                   | Contents                                  |
|----------------------------------------|-------------------------------------------|
| `<channel>_analysis_with_zscore_.csv`  | All cells with features and z-scores      |
| `<channel>_analysis_in_range_.csv`     | In-range cells only                       |
| `<channel>_analysis_out_of_range_.csv` | Out-of-range cells only                   |

In the pooled outputs, the `batch_id` column (for example `Q23_0`) identifies
the experiment each cell came from, and `channel_id` identifies the channel.

## Repository contents

This repository holds two independent analyses of huntingtin and Golgi
structure. They share subject matter but not data, design, or code.

| Path | Analysis |
|------|----------|
| `analyse_sub_folders.py` | Post-processing and outlier screening of per-cell Golgi morphology measurements (genotypes `KO`, `HAP40KO`, `Q23`, `Q145`). Described above. |
| `scatter-analysis/` | Dual-channel protein scatter analysis of fluorescent microscopy images across a 2-bromopalmitate treatment and washout timecourse. See [`scatter-analysis/README.md`](scatter-analysis/README.md). |

Note that "channel" means different things in the two analyses: an
experimental genotype group above, and a fluorescence emission channel
(C1/C2/C3) in `scatter-analysis/`.
