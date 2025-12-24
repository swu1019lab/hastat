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

### Requirements
- Python 3.9+
- Dependencies: `pandas`, `numpy`, `scipy`, `statsmodels`, `scikit-allel`, `gffutils`, `pysam`, `matplotlib`, `tomli`, `prettytable`, `networkx`

### Install via Source
```bash
git clone https://github.com/swu1019lab/hastat.git
cd hastat

# Option 1: Using build (Recommended)
pip install build --user
python -m build
pip install dist/hastat-1.0.0.tar.gz --user

# Option 2: Using setup.py
python setup.py install --user
```

## 🚀 Usage

The general syntax for `hastat` is:
```bash
hastat <subcommand> [options]
```

### Quick Start Examples

#### 1. View Haplotypes
Extract haplotype groups for a specific gene:
```bash
hastat view -v data.vcf.gz -a annotation.gff3 -i GeneID -t group -o output_prefix
```

#### 2. Statistical Analysis
Test if haplotypes are associated with phenotypic traits:
```bash
hastat stat -g hap_groups.csv -p phenotypes.csv -o stats_results.csv
```

#### 3. Visualization
Plot haplotype frequencies as a bar chart:
```bash
hastat plot bar --sample_hap hap_groups.csv --sample_group pop_info.csv -o plot.png
```

---

## 📚 Documentation

<details>
<summary><strong>1. View Module (`hastat view`)</strong></summary>

Analyze haplotypes and genetic statistics.

| Argument | Description |
|:---|:---|
| `-v, --vcf` | Input VCF file(s). Required. |
| `-a, --gff` | GFF3 annotation file. |
| `-i, --gene_id` | Target Gene ID. |
| `-r, --region` | Target region (`chr:start-end`). |
| `--homo` | Comma-separated homologous gene IDs. |
| `-t, --type` | Analysis type: `geno`, `table`, `group`, `freq`, `pi`, `fst`, `compare`. |
| `-g, --group` | Population group CSV file. |
| `-w, --size` | Window size for sliding window analysis (default: 1). |
| `-s, --step` | Step size for sliding window analysis (default: 1). |
| `--het` | Heterozygosity filter threshold (0-1). |
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

## 💡 Advanced Examples

### 🧬 Multi-Population Comparison
Compare haplotype frequencies and diversity between populations.

```bash
# 1. Generate comparison data
hastat view -t compare -v pop1.vcf.gz pop2.vcf.gz -n Pop1 Pop2 -a ann.gff -i GeneID -o comparison

# 2. Statistical analysis
hastat stat -g comparison.merged.group.csv -p traits.csv -o stat_results.csv

# 3. Visualize
hastat plot bar --sample_hap comparison.merged.group.csv --sample_group pop_info.csv -o comparison_bar.png
```

### 🕸️ Haplotype Network Analysis
Create a haplotype network to visualize evolutionary relationships.

```bash
# 1. Generate haplotype table
hastat view -v data.vcf -a ann.gff -i GeneID -t table -o hap_table

# 2. Construct network (MST)
hastat network -t hap_table.csv -o network.txt -m MST

# 3. Plot network
hastat plot network \
    --file network.txt \
    --sample_hap hap_groups.csv \
    --sample_group pop_info.csv \
    --layout spring \
    --node_color '#C5504B' '#FCE988' '#90CAEE' \
    -o network_plot.png
```

### 📉 Selection Sweep Analysis ($\pi$ & $F_{ST}$)
Scan for selection signatures using sliding windows.

```bash
# Nucleotide Diversity (Pi)
hastat view -t pi -v data.vcf -a ann.gff -i GeneID \
    -g pop_info.csv --size 1000 --step 100 -o pi_results

# Fixation Index (Fst)
hastat view -t fst -v data.vcf -a ann.gff -i GeneID \
    -g pop_info.csv --size 1000 --step 100 -o fst_results
```

## 📄 Citation

If you use **hastat** in your research, please cite:

> **Xiaodong Li, & Kun Lu. (2024).** hastat: a python library to perform gene haplotypes analysis in natural populations. *Zenodo*. https://doi.org/10.5281/zenodo.11183815