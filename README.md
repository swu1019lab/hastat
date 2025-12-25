<div align="center">

<img src="tests/hastat.png" alt="hastat logo" width="200"/>

**A Comprehensive Toolkit for Gene Haplotype Analysis in Natural Populations**

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.11001623-blue)](https://doi.org/10.5281/zenodo.11183815)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-BSD3-green)](LICENSE)

[Features](#-features) • [Installation](#-installation) • [Usage](#-usage) • [Documentation](#-documentation) • [Citation](#-citation)

</div>

---

## 📖 Description

**hastat** is a powerful Python library designed for geneticists and evolutionary biologists. It streamlines the process of viewing, analyzing, and visualizing genetic variation data, specifically focusing on gene haplotypes in natural populations. From raw VCF files to publication-ready visualizations, `hastat` provides a complete workflow.

## ✨ Features

### Core Modules
| Module | Description |
|:---|:---|
| 🧬 **`view`** | Extract and analyze genotype, haplotype, $\pi$, $F_{ST}$, and LD data from VCF files. |
| 📊 **`stat`** | Perform robust statistical analysis (ANOVA, multiple comparisons) linking haplotypes to phenotypes. |
| 🎨 **`plot`** | Generate high-quality visualizations: Bar, Pie, Box, Network, and Gene haplotype plots. |
| 🕸️ **`network`** | Construct haplotype networks using Minimum Spanning Tree (MST) and Minimum Spanning Network (MSN). |
| 🔬 **`gwas`** | Integrated wrapper for GEMMA to perform Genome-Wide Association Studies. |

### Advanced Capabilities
- ✅ **Multi-population Analysis**: Compare haplotype frequencies and diversity across different groups.
- ✅ **Selection Sweeps**: Detect signatures of selection using nucleotide diversity ($\pi$) and fixation index ($F_{ST}$).
- ✅ **Homologous Genes**: Analyze multiple homologous genes simultaneously.
- ✅ **Customization**: Flexible grouping, filtering (heterozygosity), and statistical methods.

## 📦 Installation

### 1. Prerequisites
- **Python 3.9+**
- **Dependencies**: `pandas`, `numpy`, `scipy`, `statsmodels`, `scikit-allel`, `gffutils`, `pysam`, `matplotlib`, `tomli`, `prettytable`, `networkx`

You can check your Python version with:
```bash
python --version
```

### 2. Download Source Code
Clone the repository from GitHub to your local machine. **Important: You must enter the `hastat` directory before installation.**
```bash
git clone https://github.com/swu1019lab/hastat.git
cd hastat
```

### 3. Install hastat
You can install `hastat` using either the modern `build` system (recommended) or `setup.py`.

**Tip for users in China:** To speed up dependency downloading, you can use the Tsinghua PyPI mirror by adding `-i https://pypi.tuna.tsinghua.edu.cn/simple` to pip commands.

#### Option A: Using `build` (Recommended)
This method builds a standard distribution package and installs it.

1.  **Install the build tool** (skip if already installed):
    ```bash
    pip install build --user -i https://pypi.tuna.tsinghua.edu.cn/simple
    ```
2.  **Build the package**:
    ```bash
    python -m build
    ```
    *Explanation: This command compiles the source code and creates a `dist/` directory containing the installable package file (e.g., `hastat-1.0.0.tar.gz`).*
3.  **Install the package**:
    ```bash
    pip install dist/hastat-1.0.0.tar.gz --user -i https://pypi.tuna.tsinghua.edu.cn/simple
    ```
    *Explanation: This installs the `hastat` library and its dependencies (pandas, numpy, etc.) into your Python environment.*

#### Option B: Using `setup.py` (Legacy)
```bash
python setup.py install --user
```

### 4. Installation Location
When using the `--user` flag (recommended to avoid permission issues):
- **Linux/macOS**: The executable `hastat` is typically placed in `~/.local/bin`.
  - *Note: You may need to add this path to your `$PATH` environment variable if `hastat` command is not found.*
- **Windows**: It is placed in your user's Python Scripts directory, e.g., `C:\Users\YourName\AppData\Roaming\Python\Python39\Scripts`.

To verify the installation, run:
```bash
hastat --help
```

## 🚀 Usage

The general syntax for `hastat` is:
```bash
hastat <subcommand> [options]
```

Here is a step-by-step guide based on a typical workflow:

### Step 1: Generate Haplotype Table
Extract haplotypes for a specific gene from a VCF file.
```bash
hastat view \
    -t table \
    -v data.vcf.gz \
    -a annotation.gff3 \
    -i GeneID \
    -u 2000 \
    -o results/GeneID.table
```

### Step 2: Group Samples by Haplotype
Group samples based on their haplotypes.
```bash
hastat view \
    -t group \
    -v data.vcf.gz \
    -a annotation.gff3 \
    -i GeneID \
    -u 2000 \
    -o results/GeneID.group
```

### Step 3: Statistical Analysis
Perform statistical tests to associate haplotypes with phenotypes.
```bash
hastat stat \
    -g results/GeneID.group.csv \
    -p phenotype.csv \
    -a GeneID \
    -o results/GeneID.stat.csv
```

### Step 4: Visualization

#### Haplotype Frequency (Pie Chart)
```bash
hastat plot pie \
    --sample_hap results/GeneID.group.csv \
    --sample_group sample_groups.csv \
    -o results/GeneID.pie.pdf
```

#### Phenotype Distribution (Box Plot)
```bash
hastat plot box \
    --sample_hap results/GeneID.group.csv \
    --sample_phe phenotype.csv \
    --phe_index 1 \
    -o results/GeneID.box.pdf
```

#### Haplotype Network
```bash
# 1. Generate network file
hastat network -t results/GeneID.table.csv -o results/GeneID.network.txt

# 2. Plot network
hastat plot network \
    --file results/GeneID.network.txt \
    --sample_hap results/GeneID.group.csv \
    --sample_group sample_groups.csv \
    -o results/GeneID.network.pdf
```

#### Gene Structure and Selection Signals
Visualize gene structure along with $\pi$ and $F_{ST}$ tracks.
```bash
hastat plot gene \
    --gff annotation.gff3 \
    --genes GeneID \
    --toml gene.toml \
    -o results/GeneID.gene.pdf \
    --upstream 2000
```

### Step 5: Selection Sweep Analysis
Calculate $F_{ST}$ and $\pi$ using a sliding window approach.
```bash
# Calculate Fst
hastat view \
    -t fst \
    -v data.vcf.gz \
    -a annotation.gff3 \
    -i GeneID \
    -u 1000000 \
    -d 1000000 \
    -g sample_groups.csv \
    -o results/GeneID.fst \
    --size 10000 \
    --step 1000

# Calculate Pi
hastat view \
    -t pi \
    -v data.vcf.gz \
    -a annotation.gff3 \
    -i GeneID \
    -u 1000000 \
    -d 1000000 \
    -g sample_groups.csv \
    -o results/GeneID.pi \
    --size 10000 \
    --step 1000
```

### 1. View Module (`hastat view`)
The `view` module is the core tool for extracting genetic data and calculating population genetics statistics.

#### 📍 Flexible Input Modes
`hastat` allows you to define target loci in multiple ways:
- **Single Gene**: `-i GeneID` (Requires GFF file via `-a`)
- **Genomic Region**: `-r chr:start-end` (e.g., `chr1:1000-5000`)
- **Batch List**: `-l gene_list.txt` (Process a list of gene IDs automatically)
- **Homologous Genes**: `--homo GeneA,GeneB` (Analyze multiple homologs as a single unit)

#### 📏 Region Extension & Control
You can extend the analysis scope beyond the gene body, which is crucial for analyzing promoter regions or regulatory elements.
- **`--upstream <int>`**: Extend $N$ bp upstream of the gene/region start.
- **`--downstream <int>`**: Extend $N$ bp downstream of the gene/region end.

> **Note**: These parameters apply to **all** analysis types. For example, if you calculate $F_{ST}$ with `-u 2000`, the sliding window analysis will start 2kb upstream of the gene.

#### 🧬 Analysis Types (`-t/--type`)
| Type | Description |
|:---|:---|
| `geno` | Extract raw genotype matrix (0/0, 0/1, 1/1). |
| `table` | Generate a haplotype table for samples. |
| `group` | Group samples by their haplotypes. |
| `freq` | Calculate haplotype frequencies. |
| `compare` | **Multi-Population Mode**: Compare haplotypes across multiple VCF files. |
| `pi` | Calculate Nucleotide Diversity ($\pi$). |
| `fst` | Calculate Fixation Index ($F_{ST}$). |

#### 📊 Population Genetics ($\pi$ & $F_{ST}$)
`hastat` calculates $\pi$ and $F_{ST}$ using a **sliding window** approach.

- **Method**:
    1.  **Region Definition**: The total analysis region is determined by the gene coordinates plus any `--upstream` or `--downstream` extension.
    2.  **Sliding Windows**: Windows of size `--size` move across this region with a step of `--step`.
    3.  **Calculation**: Metrics are computed for each window based on sample groups provided in `-g`.
- **Required Arguments**:
    - `-g, --group`: CSV file defining sample populations (Columns: Sample, Group).
    - `-w, --size`: Window size in bp (e.g., `1000`).
    - `-s, --step`: Step size in bp (e.g., `100`).

**Example: Selection Sweep Analysis**
Calculate $F_{ST}$ for a gene and its surrounding region (extended 1MB upstream and downstream) using a sliding window approach (Reference: [Nature Genetics, 2020](https://www.nature.com/articles/s41588-020-0604-7)):
```bash
hastat view -v data.vcf.gz -a ann.gff3 -i GeneID \
    -t fst -g sample_groups.csv \
    -u 1000000 -d 1000000 --size 10000 --step 10000 -o output_prefix
```

### 2. Statistics Module (`hastat stat`)
Perform statistical tests to associate haplotypes with phenotypes.
- **Input**: Haplotype groups (from `view -t group`) and Phenotype data (`-p`).
- **Method**: ANOVA followed by multiple comparisons (`TukeyHSD` or `AllPairTest`).

```bash
hastat stat -g hap_groups.csv -p traits.csv -a GeneID -o stat_results.csv
```

### 3. Visualization (`hastat plot`)
Generate high-quality visualizations for your data.

#### Common Plot Types
- **`bar` / `pie`**: Visualize haplotype frequencies in different populations.
- **`box`**: Compare phenotypic values among haplotypes or populations.
    - Supports significance testing (t-test, U-test).
    - Use `--comparisons` to specify pairs for testing.

#### Gene Structure & Tracks (`hastat plot gene`)
Visualize the gene model along with haplotype variations and statistical tracks ($\pi$, $F_{ST}$).
- **Configuration**: Uses a TOML file to control layout, colors, and tracks.
- **Features**:
    - Plot CDS/Exon structure.
    - Align haplotype variations to the gene.
    - Display sliding window tracks for selection signatures.

**Example: Plotting Gene with Tracks**
```bash
hastat plot gene --gff ann.gff3 --genes GeneID --toml config.toml \
    --upstream 2000 --downstream 1000 -o gene_plot.pdf
```

### 4. Network Analysis (`hastat network`)
Construct and visualize haplotype networks (Minimum Spanning Tree/Network).
1.  **Generate Table**: `hastat view -t table ...`
2.  **Build Network**: `hastat network -t table.csv -o network.txt`
3.  **Visualize**: `hastat plot network --file network.txt ...`

## 💡 Advanced Examples

### Multi-Population Comparison
Compare haplotype frequencies between two populations (e.g., Panel1 vs Panel2).

```bash
# 1. Compare haplotypes from two VCFs
hastat view -t compare \
    -v panel1.vcf.gz panel2.vcf.gz \
    -n Panel1 Panel2 \
    -i GeneID -a ann.gff3 -o comparison_result

# 2. Visualize the comparison
hastat plot bar --sample_hap comparison_result.merged.group.csv \
    --sample_group pop_info.csv -o comparison_bar.pdf
```

### Large-Scale Selection Sweep
Scan a large genomic region (e.g., 1Mb upstream/downstream) for selection signals.

```bash
hastat view -t fst \
    -v data.vcf.gz -a ann.gff3 -i GeneID \
    -u 1000000 -d 1000000 \
    -g sample_groups.csv \
    --size 100000 --step 10000 \
    -o large_scale_fst
```

## 📚 Documentation

<details>
<summary><strong>1. View Module (`hastat view`)</strong></summary>

The `view` module is the entry point for extracting and analyzing genetic data. It supports flexible input methods and various analysis types.

#### Input Modes
| Mode | Argument | Description |
|:---|:---|:---|
| **Single Gene** | `-i GeneID` | Analyze a specific gene. Requires `-a/--gff`. |
| **Region** | `-r chr:start-end` | Analyze a specific genomic region (e.g., `chr1:1000-2000`). |
| **Batch Genes** | `-l list.txt` | Analyze a list of genes provided in a file (one ID per line). Requires `-a/--gff`. |
| **Homologous** | `--homo ID1,ID2` | Analyze multiple homologous genes simultaneously. Requires `-a/--gff`. |

#### Region Control
You can extend the analysis region (e.g., to include promoters) using:
- `-u, --upstream <int>`: Extend $N$ bp upstream.
- `-d, --downstream <int>`: Extend $N$ bp downstream.

#### Analysis Types (`-t, --type`)
- **`geno`**: Extract genotype matrix.
- **`table`**: Generate haplotype table.
- **`group`**: Group samples by haplotypes.
- **`freq`**: Calculate haplotype frequencies.
- **`compare`**: Compare haplotypes across multiple populations (VCFs).
- **`pi`**: Calculate Nucleotide Diversity ($\pi$).
- **`fst`**: Calculate Fixation Index ($F_{ST}$).

#### Population Genetics Statistics ($\pi$ & $F_{ST}$)
For `pi` and `fst` analysis, `hastat` uses a sliding window approach:
- **Windowing**: Controlled by `-w/--size` (window size, default 1bp) and `-s/--step` (step size, default 1bp).
- **Grouping**: Requires a sample group file via `-g/--group` (csv: sample,group).
- **Calculation**:
    - **$\pi$**: Calculated for each group in the window.
    - **$F_{ST}$**: Calculated for all pairwise combinations of groups in the window.
    - **Region**: The analysis region is determined by the gene/region input plus any upstream/downstream extension.

**Example:** Calculate $F_{ST}$ for a gene plus 2kb upstream, using 1kb windows with 100bp step:
```bash
hastat view -v data.vcf.gz -a ann.gff3 -i GeneID -u 2000 -t fst -g groups.csv --size 1000 --step 100
```
</details>

<details>
<summary><strong>2. Statistics Module (`hastat stat`)</strong></summary>

Perform statistical tests (ANOVA, TukeyHSD).

| Argument | Description |
|:---|:---|
| `-g, --group` | Haplotype group CSV file. Required. |
| `-p, --pheno` | Phenotype CSV file. Required. |
| `-m, --method` | Method: `TukeyHSD` (default) or `AllPairTest`. |
| `-s, --min_hap_size` | Minimum sample size per haplotype (default: 10). |
</details>

<details>
<summary><strong>3. Network Module (`hastat network`)</strong></summary>

Construct haplotype networks.

| Argument | Description |
|:---|:---|
| `-t, --hap_table` | Haplotype table from `view -t table`. Required. |
| `-m, --method` | `MST` (Minimum Spanning Tree) or `MSN`. |
| `-f, --threshold_factor` | Threshold factor for MSN (default: 1.2). |
</details>

<details>
<summary><strong>4. Plot Module (`hastat plot`)</strong></summary>

Visualize results. Subcommands: `bar`, `pie`, `box`, `network`, `gene`.

**Common Arguments:**
- `-o, --output`: Output file path.
- `--width`, `--height`: Figure dimensions.

**Subcommand Specifics:**
- **`bar` / `pie`**: Requires `--sample_hap`, `--sample_group`.
- **`box`**: Requires `--sample_hap`, `--sample_phe`, `--comparisons`.
- **`network`**: Requires `--file` (from network module). Supports layout customization (`--layout`, `--node_color`, etc.).
- **`gene`**: Requires `--gff`, `--genes`, `--toml`.

**Gene Plot Configuration (`gene.toml`):**
Control the appearance of gene structure plots using a TOML file.

| Parameter | Type | Description | Default |
|:---|:---|:---|:---|
| **Gene Structure** | | | |
| `bbox` | list | Bounding box `[x, y, width, height]` relative to figure. | `[0.1, 0.1, 0.8, 0.8]` |
| `xy` | list | Gene anchor point `[x, y]` relative to bbox. | `[0, 0.5]` |
| `width` | float | Gene width relative to bbox. | `1.0` |
| `height` | float | Gene height relative to bbox. | `0.1` |
| `color` | string | Gene color. | `"C1"` |
| `feature` | string | Feature to plot (e.g., "CDS", "exon"). | `"CDS"` |
| `gene_label_show` | bool | Whether to show gene name label. | `true` |
| `gene_label_size` | int | Font size for gene label. | `10` |
| **Annotation** | | | |
| `ann_file` | string | Path to variant annotation CSV file. | `""` |
| `ann_y` | float | Y position for annotations relative to bbox. | `0.8` |
| `ann_text_size` | int | Font size for annotation text. | `8` |
| **Haplotypes** | | | |
| `hap_file` | string | Path to haplotype data CSV file. | `""` |
| `hap_y` | float | Y position for haplotype matrix relative to bbox. | `0.0` |
| `hap_height` | float | Height of haplotype matrix relative to bbox. | `0.3` |
| `hap_cols` | list | List of specific haplotypes to show (empty for all). | `[]` |
| `hap_label_show` | bool | Whether to show haplotype labels. | `true` |
| `hap_label_size` | int | Font size for haplotype labels. | `8` |
| **Haplotype Groups** | | | |
| `hap_group_file` | string | Path to haplotype group CSV file. | `""` |
| `hap_group_bar_x` | float | X position offset for group bar. | `-0.12` |
| `hap_group_bar_width` | float | Width of group bar. | `0.02` |
| `hap_freq_show` | bool | Whether to show frequency tracks. | `true` |
| `hap_freq_height` | float | Height of each frequency track. | `0.05` |
| `hap_freq_y` | float | Y position for frequency tracks. | `-0.2` |
| `hap_freq_plot_type` | string | Plot type: `"bar"` or `"pie"`. | `"bar"` |
| **Tracks (Pi/Fst)** | | | |
| `pi_file` | string | Path to Pi data CSV file. | `""` |
| `pi_y` | float | Y position for Pi plot. | `0.6` |
| `pi_height` | float | Height of Pi plot. | `0.15` |
| `pi_style` | string | Plot style: `"line"` or `"fill"`. | `"fill"` |
| `fst_file` | string | Path to Fst data CSV file. | `""` |
| `fst_y` | float | Y position for Fst plot. | `0.8` |
| `fst_height` | float | Height of Fst plot. | `0.15` |
| `fst_style` | string | Plot style: `"line"` or `"fill"`. | `"fill"` |

```toml
[default]
# Example configuration
bbox = [0.1, 0.1, 0.8, 0.8]
color = "C1"
gene_label_show = true
hap_file = "haplotypes.csv"
pi_file = "pi.csv"
```
</details>

<details>
<summary><strong>5. GWAS Module (`hastat gwas`)</strong></summary>

Run GWAS using GEMMA.

| Argument | Description |
|:---|:---|
| `-c, --config` | Configuration TOML file. Required. |
</details>

---

## 📄 Citation

If you use **hastat** in your research, please cite:

> **Xiaodong Li, & Kun Lu. (2024).** hastat: a python library to perform gene haplotypes analysis in natural populations. *Zenodo*. https://doi.org/10.5281/zenodo.11183815
