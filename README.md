<div align="center">

<img src="https://raw.githubusercontent.com/swu1019lab/hastat/master/docs/assets/hastat.png" alt="hastat logo" width="200"/>

**A Comprehensive Toolkit for Gene Haplotype Analysis in Natural Populations**

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.11001623-blue)](https://doi.org/10.5281/zenodo.11183815)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-BSD3-green)](LICENSE)

[Features](#-features) • [Installation](#-installation) • [Example results](#-example-results) • [Usage](#-usage) • [Documentation](#-documentation) • [Citation](#-citation)

</div>

---

## 📖 Description

**hastat** is a powerful Python library designed for geneticists and evolutionary biologists. It streamlines the process of viewing, analyzing, and visualizing genetic variation data, specifically focusing on gene haplotypes in natural populations. From raw VCF files to publication-ready visualizations, `hastat` provides a complete workflow.

## ✨ Features

### Core Modules
| Module | Description |
|:---|:---|
| 🧬 **`view`** | Extract and analyze genotype, haplotype, $\pi$, $F_{ST}$ and Tajima's D data from VCF files. LD is available through the Python API. |
| 📊 **`stat`** | Perform robust statistical analysis (ANOVA, multiple comparisons) linking haplotypes to phenotypes. |
| 🎨 **`plot`** | Generate high-quality visualizations: Bar, Pie, Box, Network, and Gene haplotype plots. |
| 🕸️ **`network`** | Construct haplotype networks using Minimum Spanning Tree (MST) and Minimum Spanning Network (MSN). |

### Advanced Capabilities
- ✅ **Multi-population Analysis**: Compare haplotype frequencies and diversity across different groups.
- ✅ **Selection Signals Scan**: Detect signatures of selection using nucleotide diversity ($\pi$), fixation index ($F_{ST}$), and Tajima's D.
- ✅ **Homologous Genes Analysis**: Analyze multiple homologous genes simultaneously.
- ✅ **Customization parameters**: Flexible grouping, filtering (heterozygosity), and statistical methods.

### Comparison with Other Tools
| Feature | hastat | plink | vcftools | bcftools |
| :--- | :--- | :--- | :--- | :--- |
| **Tool Positioning** | Genotype Manipulation and Gene-centric Haplotype Analysis | Genotype Manipulation and GWAS | Genotype Manipulation and Statistics | Genotype Manipulation and Variant Calling |
| **Major Functions** | Mining Favorable Haplotypes and Breeding Improvement | Association testing, LD, Filtering | Filtering, Calc Pi/Fst | Filtering, Merging, Calling |
| **Analysis Range** | Gene and nearby regions | Whole genome | Whole genome | Whole genome |
| **Primary Input** | Genotype | Genotype | Genotype | Genotype |
| **Main Output** | Haplotype | Genotype | Genotype | Genotype |
| **Genotype fetch** | YES | YES | YES | YES |
| **Haplotype fetch** | YES | NO | NO | NO |
| **Phenotype Association** | YES | YES | NO | NO |
| **Haplotype Visuals** | YES | NO | NO | NO |
| **Multiple Species** | YES | YES | YES | YES |


### ⚠️ Limitations
`hastat` is specialized for gene-centric haplotype analysis and visualization. It is **NOT** designed for:
- **Raw Data Processing**: It does not handle FASTQ/BAM files or perform variant calling.
- **Whole-Genome Operations**: For filtering or merging VCFs, use `bcftools` or `vcftools`.
- **Complex Pedigrees**: It focuses on natural populations, not family-based linkage analysis.

## 📦 Installation

### 1. Prerequisites
> **⚠️ Windows Compatibility**: Since `pysam` is not available for Windows, `hastat` **cannot** be installed directly on Windows. Windows users are highly recommended to use **WSL (Windows Subsystem for Linux)**.

> **✅ Tested Environments**: We have tested the installation and download on **CentOS Linux 8** and **Ubuntu 24**.

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
You can install `hastat` using the modern `build` system or directly with pip.

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
    *This creates a `dist/` directory containing the source archive and wheel.*
3.  **Install the package**:
    ```bash
    python -m pip install dist/hastat-2.0.0-py3-none-any.whl --user -i https://pypi.tuna.tsinghua.edu.cn/simple
    ```
    *Explanation: This installs the `hastat` library and its dependencies (pandas, numpy, etc.) into your Python environment.*

#### Option B: Install this source checkout directly
```bash
python -m pip install . --user
```

### Input table conventions

Sample and group CSVs require a header and one row per unique sample. Identifiers
and group labels are read as text: `001`, `1` and `NA` remain distinct. Surrounding
whitespace is stripped; empty identifiers/labels and duplicate IDs are rejected.
Population assignment files have exactly two columns, and every listed sample
must exist in the VCF. Phenotype tables use a `samples` column plus numeric trait
columns; missing tokens such as `NA` apply only to trait values. Network sample
joins use matched samples and warn about unmatched rows; duplicate assignments
and group indexes outside the available columns fail explicitly.

### 4. Installation Location
When using the `--user` flag (recommended to avoid permission issues):
- **Linux/macOS**: The executable `hastat` is typically placed in `~/.local/bin`.
  - *Note: You may need to add this path to your `$PATH` environment variable if `hastat` command is not found.*

To verify the installation, run:
```bash
hastat --help
```

## 🖼️ Example results

Examples use 160 synthetic samples and eight haplotypes. See the
[gallery](https://swu1019lab.github.io/hastat/gallery.html) and [example guide](https://swu1019lab.github.io/hastat/examples.html) for inputs and commands.

### Gene haplotypes with 1 Mb flanks

Gene structure, haplotypes, variant annotations, allele frequencies and LD,
with π, Fst and Tajima’s D tracks spanning 1 Mb flanks.

[Recommended TOML](docs/examples/gene_tracks.toml) · [Reproduction commands](docs/examples/commands.sh)

![Gene and haplotypes aligned with Pi, Fst and Tajima's D windows](docs/images/gene_tracks.png)

### MST and MSN haplotype networks

MST and MSN networks showing sample counts, population composition and differences between haplotypes.

![Minimum spanning tree and minimum spanning network of the eight haplotypes](docs/images/mst_msn.png)

## 🚀 Usage

The general syntax for `hastat` is:
```bash
hastat <subcommand> [options]
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
- **`--up <int>`**: Extend $N$ bp upstream of the gene/region start.
- **`--down <int>`**: Extend $N$ bp downstream of the gene/region end.

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
| `tajima_d` | Calculate Tajima's D. |

For single-VCF `table`, `group` and `freq`, all samples use the same complete
loci after the requested genotype filters. The three outputs share one
deterministic profile/rank mapping, including multi-region requests. `freq`
contains integer sample counts: `haplotypes,samples` without metadata, or
`haplotypes,<population names>` with a single header when grouped. Group files
may select a sample subset for counting; they do not redefine the profiles.
Duplicate or unknown sample assignments are rejected. Empty outputs include column headers.

#### 📊 Population Genetics ($\pi$, $F_{ST}$ & Tajima's D)

`hastat` calculates $\pi$, $F_{ST}$ and Tajima's D using a **sliding window** approach.

- **Nucleotide Diversity ($\pi$)**: Measures the average number of nucleotide differences per site between two randomly chosen sequences in a population. It reflects the genetic diversity within a population. High $\pi$ indicates high genetic diversity.
- **Fixation Index ($F_{ST}$)**: Estimates differentiation between the selected populations using Weir–Cockerham variance components. Finite-sample estimates can be negative and are retained; undefined estimates are NaN. Window estimates use a ratio of component sums, not a mean of per-SNP Fst values.
- **Tajima's D**: A test of neutrality. It compares two measures of nucleotide diversity ($\pi$ and $\theta_W$).
    - $D \approx 0$: Neutral evolution.
    - $D > 0$: Balancing selection or sudden population contraction.
    - $D < 0$: Positive selection (selective sweep) or population expansion.

- **Method**:
    1.  **Region Definition**: The total analysis region is determined by the gene coordinates plus any `--up` or `--down` extension.
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

Fst uses the same 1-based inclusive window boundaries as pi/Tajima's D. Group
assignments must refer to existing, unique samples and contain at least two
populations. Records lacking enough called individuals for the variance
components are excluded; empty or unestimable windows retain NaN. For pi/Fst/Tajima's D, use the default `--maf` and `--max_missing` settings.
The Fst CSV format is `chrom,start,end,counts,fst,groups`;
`counts` is the input variant count.

### 2. Statistics Module (`hastat stat`)
Perform statistical tests to associate haplotypes with phenotypes.
- **Input**: Haplotype groups (from `view -t group`) and Phenotype data (`-p`).
- **Methods**: `TukeyHSD` (default): classical equal-variance ANOVA plus Tukey–Kramer comparisons; `WelchBonferroni`: Welch ANOVA plus two-sided Welch t tests with Bonferroni correction.

#### Phenotype Data Format
The phenotype data file should be a CSV file with the following format:
- **First column**: Sample ID (must match sample names in haplotype group file).
- **Remaining columns**: Phenotypic traits. You can have multiple traits (columns), and each will be analyzed separately.
- **Header**: First row should contain column names (samples, trait1, trait2, ...).

**Example phenotype data (`phenotype.csv`):**
```
samples,Height,Weight,Yield,Quality
Sample_1,175.5,65.2,8500.5,7.8
Sample_2,172.3,62.8,8200.3,7.5
Sample_3,178.9,68.1,9100.2,8.1
Sample_4,170.5,61.5,7800.8,7.2
Sample_5,176.2,66.5,8900.5,8.0
```

**Usage:**
```bash
hastat stat -g hap_groups.csv -p phenotype.csv -a GeneID -o stat_results.csv
```

For Welch ANOVA and Bonferroni-corrected Welch comparisons:
```bash
hastat stat -g hap_groups.csv -p phenotype.csv -m WelchBonferroni -o stat_results.csv
```

**Tool outputs:** `-o stat_results.csv` writes only the requested full-precision
pairwise CSV table. Its fields include the method, sample counts, group means
and standard deviations, mean differences, confidence intervals, adjusted P
values, rejection decisions, omnibus P values, and comparison-family sizes.
Without `-o`, the tool prints the comparison table and omnibus records to stdout.

### 3. Visualization (`hastat plot`)
Generate high-quality visualizations for your data.

The Python LD APIs `hastat.utils.population.calculate_ld` / `calculate_multi_loci_ld`
take a `VariantData` object as their first argument and return
square r² matrices: valid diagonals are 1, undefined entries are NaN and empty
regions return (0,0). Multi-region overlaps count once. LD uses pairwise complete
diploid dosage observations with the existing Rogers-Huff implementation.
`view --type ld` is not a CLI option.

```python
from hastat.dataio.vcf import read
from hastat.utils.variant import VariantData
from hastat.utils.population import calculate_ld

with read("data.vcf.gz") as vcf:
    variants = VariantData(vcf)
    r_squared = calculate_ld(variants, "chr1", 1000, 5000)
```

For haplotype extraction, use `hastat.utils.haplotype.extract_haplotypes(variants, loci)`
to obtain the haplotype table and sample groups. Haplotype–phenotype association is
available through `hastat.utils.association.HaplotypeAssociation`:

```python
from hastat.dataio.table import read_sample_csv
from hastat.utils.association import HaplotypeAssociation

phenotypes = read_sample_csv("phenotype.csv")
groups = read_sample_csv("hap_groups.csv")
analysis = HaplotypeAssociation(phenotypes, groups, method="WelchBonferroni")
comparisons = analysis.comparisons
```

#### Common Plot Types
- **`bar` / `pie`**: Visualize haplotype frequencies in different populations.
- **`box`**: Compare phenotypic values among haplotypes or populations.
    - `TukeyHSD` and `WelchBonferroni` reuse the `stat` implementation.
      The existing default `t-test` and `mannwhitneyu` use Bonferroni-adjusted
      two-sided P values; `welch` uses `WelchBonferroni`.
    - `--comparisons` selects brackets, not the multiple-testing family.
      The family includes all eligible groups in the input, even when
      `--haplotypes` displays a subset. `--min_hap_size` defaults to 2 for plots
      (10 for `stat`); use identical inputs, method and minimum for agreement.
    - Undefined tests have no bracket; the reason is logged. `ns` is reserved
      for finite adjusted P values >= 0.05. No extra output table is created.

Box/bar/pie sample identifiers and group labels are read as text, preserving
`001`, `1`, and `NA`. Duplicate/empty IDs and duplicate headers are errors.
Input tables are joined one-to-one on shared IDs; unmatched rows are logged.
Invalid plotting input produces a nonzero CLI exit.

#### Gene Structure & Tracks (`hastat plot gene`)
Visualize the gene model along with haplotype variations and statistical tracks ($\pi$, $F_{ST}$).
- **Configuration**: Uses a TOML file to control layout, colors, and tracks.
- **Features**:
    - Plot CDS/Exon structure.
    - Align haplotype variations to the gene.
    - Display sliding window tracks for selection signatures.

**Example: Plotting Gene with Tracks**
```bash
hastat plot gene --gff ann.gff3 --gene GeneID --toml config.toml \
    --up 2000 --down 1000 -o gene_plot.pdf
```

### 4. Network Analysis (`hastat network`)
Construct and visualize haplotype networks (Minimum Spanning Tree/Network).
1.  **Generate Table**: `hastat view -t table ...`
2.  **Build Network**: `hastat network -t table.csv -o network.txt`
3.  **Visualize**: `hastat plot network --file network.txt ...`

---

### Here is a step-by-step guide based on a typical workflow

### Step 1: Generate Haplotype Table
Extract haplotypes for a specific gene from a VCF file.
```bash
hastat view \
    -t table \
    -v data.vcf.gz \
    -a annotation.gff3 \
    -i GeneID \
    -u 2000 \
    -o results/GeneID
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
    -o results/GeneID
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
    --gene GeneID \
    --toml gene.toml \
    -o results/GeneID.gene.pdf \
    --up 2000
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
    -o results/GeneID \
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
    -o results/GeneID \
    --size 10000 \
    --step 1000
```

## 💡 Advanced Examples

### Multi-Population Comparison
Compare haplotype frequencies between two populations (e.g., Wild and Landrace).

Inputs must use the same reference and allele representation. Compare aligns
the full chromosome/position/REF/ALT key and groups samples on loci present in
every input with complete calls across all included samples. It never treats
an absent record as a reference genotype or missing calls as a haplotype state.
The existing `population_specific` field describes input-record availability,
not proof of a biologically private allele. Outputs are unphased genotype
profiles; use `MergedHap` labels for comparison across populations.

Local numbering is deterministic (abundance, then genotype tuple descending).
Sample suffixes `_1`, `_2`, ... refer to VCF input order and must match
downstream group/phenotype files. Group CSV headers are `samples,haplotypes`.
Empty outputs include column headers. Non-default `--het`, `--maf` and
`--max_missing` settings are not supported by `compare`.

```bash
# 1. Compare haplotypes from two VCFs
hastat view -t compare \
    -v wild.vcf.gz landrace.vcf.gz \
    -n Wild Landrace \
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

Read the [web documentation](https://swu1019lab.github.io/hastat/) for command-line options, input formats and gene-plot TOML settings.

Explore the [example guide](https://swu1019lab.github.io/hastat/examples.html) and [PNG gallery](https://swu1019lab.github.io/hastat/gallery.html) for complete workflows and figures.

---

## 📄 Citation

If you use **hastat** in your research, please cite:
> **Xiaodong Li, Yuhong Wang, Boyu Meng, Zhixuan Qin, Minghao Zhang, Jiaming Song, Kun Lu**, Hastat: Mining Favorable Haplotypes of Trait-Associated Genes for Breeding Improvement in Natural Populations, *Horticulture Research*, 2026;, uhag018, https://doi.org/10.1093/hr/uhag018
