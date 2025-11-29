![hastat](tests/hastat.png)

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.11001623-blue)](https://doi.org/10.5281/zenodo.11183815)

## Table of Contents
- [Table of Contents](#table-of-contents)
  - [⭐ Description](#-description)
  - [⭐ Features](#-features)
  - [⭐ Requirements](#-requirements)
  - [⭐ Installation](#-installation)
  - [⭐ Usage](#-usage)
  - [⭐ Examples](#-examples)
    - [🏷️ 1. View gene haplotypes and genotype data](#️-1-view-gene-haplotypes-and-genotype-data)
    - [🏷️ 2. Perform haplotype statistical analysis](#️-2-perform-haplotype-statistical-analysis)
    - [🏷️ 3. Visualize haplotype data](#️-3-visualize-haplotype-data)
      - [1️⃣ Bar plot](#1️⃣-bar-plot)
      - [️2️⃣ Pie plot](#️2️⃣-pie-plot)
      - [️3️⃣ Box plot](#️3️⃣-box-plot)
      - [️4️⃣ Network plot](#️4️⃣-network-plot)
      - [️5️⃣ Gene haplotypes plot](#️5️⃣-gene-haplotypes-plot)
    - [🏷️ 4. Network analysis](#️-4-network-analysis)
    - [🏷️ 5. Multi-population analysis](#️-5-multi-population-analysis)
    - [🏷️ 6. Selection sweep analysis](#️-6-selection-sweep-analysis)
    - [🏷️ 7. GWAS analysis](#️-7-gwas-analysis)
  - [⭐ Advanced Usage](#-advanced-usage)
  - [⭐ Citation](#-citation)

### ⭐ Description
A comprehensive Python library for gene haplotype analysis in natural populations, providing tools for viewing, analyzing, and visualizing genetic variation data.

### ⭐ Features
Main modules:
- [x] `view` - View genotype, haplotype, pi, fst, ld data of interested genes in VCF files
- [x] `stat` - Perform haplotype statistical analysis with multiple phenotypes using ANOVA and multiple comparisons
- [x] `plot` - Visualize haplotype data using bar, pie, box, network, and gene haplotye plots
- [x] `network` - Perform haplotype network analysis using MST and MSN methods
- [x] `gwas` - Perform GWAS analysis using GEMMA wrapper

Advanced features:
- [x] Multi-population comparison analysis
- [x] Homologous gene analysis
- [x] Selection sweep analysis (π and FST)
- [x] Custom population grouping
- [x] Heterozygosity filtering
- [x] Multiple statistical test methods

### ⭐ Requirements
Python 3.9 or higher and the following packages are required:
- pandas >= 1.3.3
- numpy >= 1.21.2
- scipy >= 1.7.1
- statsmodels >= 0.13.0
- scikit-allel >= 1.3.6
- gffutils >= 0.10.1
- pysam >= 0.17.0
- matplotlib
- tomli
- prettytable
- networkx

### ⭐ Installation

```bash
git clone https://github.com/swu1019lab/hastat.git
cd hastat
# install the package using build commands (recommended)
pip install build --user
python -m build
pip install dist/hastat-1.0.0.tar.gz --user
# or
# install the package using setup.py install command
python setup.py install --user
```

### ⭐ Usage
```bash
usage: hastat [-h] [--version] [--log LOG] {view,stat,plot,network,gwas} ...

A package for gene haplotype analysis in natural populations

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --log LOG             The log file name (default: stdout)

subcommands:
  valid subcommands

  {view,stat,plot,network,gwas}
                        additional help
    view                View and analyze the haplotypes data of genes or target regions
    stat                Perform gene statistical analysis for related-traits
    plot                Visualize gene haplotypes analysis results
    network             Perform gene haplotype network analysis
    gwas                Perform GWAS analysis using GEMMA wrapper
```

### ⭐ Examples

#### 🏷️ 1. View gene haplotypes and genotype data

**Basic haplotype group analysis:**
```bash
hastat view -v test.snps.vcf.gz -u 2000 -a test.gff3 -i gene_id -t group -o gene_id
```

**Genotype table analysis:**
```bash
hastat view -v test.vcf -a test.gff -i gene_id -t table -o gene_id
```

**Haplotype frequency analysis:**
```bash
hastat view -v test.vcf -a test.gff -i gene_id -t freq -g population_groups.csv -o gene_id
```

**Multi-population comparison:**
```bash
hastat view -t compare -v pop1.vcf.gz pop2.vcf.gz -n pop1 pop2 -a annotation.gff -i gene_id -u 2000 -o comparison_results
```

**Homologous gene analysis:**
```bash
hastat view -v test.vcf -a test.gff --homo gene1,gene2,gene3 -t group -o homologous_results
```

**Selection sweep analysis (π and FST):**

`sample_groups.csv`
```csv
samples,groups
sample1,group1
sample2,group1
sample3,group1
sample4,group1
sample5,group1
sample6,group2
sample7,group2
sample8,group2
sample9,group2
```

```bash
# π analysis
hastat view -t pi -v test.vcf -a test.gff -i gene_id -u 2000 -g sample_groups.csv -o pi_results --size 1000 --step 100

# FST analysis
hastat view -t fst -v test.vcf -a test.gff -i gene_id -u 2000 -g sample_groups.csv -o fst_results --size 1000 --step 100
```

Available view types:
- `geno` - Raw genotype data
- `table` - Haplotype table format
- `group` - Haplotype group assignments
- `freq` - Haplotype frequencies
- `pi` - Nucleotide diversity
- `fst` - Fixation index
- `compare` - Multi-population comparison

#### 🏷️ 2. Perform haplotype statistical analysis

```bash
hastat stat -g hap_groups.csv -p sample_phe.csv -o hap_stat.csv
```

**Input file formats:**

`hap_groups.csv` (haplotype assignments):
```csv
samples,haplotypes
sample1,Hap1
sample2,Hap2
sample3,Hap3
sample4,Hap4
sample5,Hap1
sample6,Hap2
sample7,Hap3
sample8,Hap4
```

`sample_phe.csv` (phenotype data):
```csv
samples,trait1,trait2,trait3
sample1,5.2,3.1,7.8
sample2,3.8,2.9,6.2
sample3,7.1,4.2,8.5
sample4,2.3,1.8,5.1
sample5,9.4,5.6,9.2
sample6,4.7,3.3,7.1
sample7,6.8,4.1,8.3
sample8,8.1,4.8,9.0
```

**Output format:**
```csv
group1,group2,meandiff,p-adj,lower,upper,reject,pheno,annotate,anova,count1,count2,mean1,mean2
Hap1,Hap2,2.0,0.7269,-3.4809,7.4809,False,trait1,gene,0.6503121273263954,5.0,5.0,4.2,6.2
Hap1,Hap3,0.8,0.9747,-4.6809,6.2809,False,trait1,gene,0.6503121273263954,5.0,5.0,4.2,5.0
```

#### 🏷️ 3. Visualize haplotype data

##### 1️⃣ Bar plot
```bash
hastat plot bar \
--sample_hap hap_groups.csv \
--sample_group sample_groups.csv \
--group_index 1 \
--calc_percentage \
--haplotypes Hap1 Hap2 Hap3 Hap4 \
--x_label "Haplotype" \
--y_label "Frequency (%)" \
-o haplotype_bar.png \
--width 10 \
--height 5
```

##### ️2️⃣ Pie plot
```bash
hastat plot pie \
--sample_hap hap_groups.csv \
--sample_group sample_groups.csv \
--group_index 1 \
--calc_percentage \
--haplotypes Hap1 Hap2 Hap3 Hap4 \
-o haplotype_pie.png \
--width 8 \
--height 6
```

##### ️3️⃣ Box plot
```bash
hastat plot box \
--sample_hap hap_groups.csv \
--sample_phe sample_phe.csv \
--phe_index 1 \
--comparisons Hap1 Hap2 \
--comparisons Hap1 Hap3 \
--comparisons Hap1 Hap4 \
--method t-test \
--haplotypes Hap1 Hap2 Hap3 Hap4 \
-o haplotype_box.png \
--width 8 \
--height 6
```

##### ️4️⃣ Network plot
```bash
# First generate network file
hastat network -t hap_table.csv -o network.txt -m MST

# Then plot the network
hastat plot network \
--file network.txt \
--sample_hap hap_groups.csv \
--sample_group sample_groups.csv \
--group_index 1 \
-o haplotype_network.png \
--show_node_label \
--node_color '#C5504B' '#FCE988' '#90CAEE' '#114F8B' \
--node_scale_factor 100 \
--layout spring \
--seed 100 \
--width 12 \
--height 8
```

##### ️5️⃣ Gene haplotypes plot
```bash
hastat plot gene \
--gff annotation.gff \
--genes gene1 \
--toml gene.toml \
--upstream 2000 \
--downstream 1000 \
-o gene_haplotypes.png \
--width 12 \
--height 8
```

#### 🏷️ 4. Network analysis

**Generate haplotype network:**
```bash
hastat network -t hap_table.csv -o network.txt -m MST
hastat network -t hap_table.csv -o network.txt -m MSN -f 1.2
```

**Network visualization:**
```bash
hastat plot network \
--file network.txt \
--sample_hap hap_groups.csv \
--sample_group sample_groups.csv \
-o network_plot.png \
--show_node_label \
--node_color '#C5504B' '#FCE988' '#90CAEE' '#114F8B' \
--node_scale_factor 100 \
--layout spring \
--seed 100
```

#### 🏷️ 5. Multi-population analysis

**Compare haplotypes across populations:**
```bash
# Generate comparison data
hastat view -t compare -v pop1.vcf.gz pop2.vcf.gz -n pop1 pop2 -a annotation.gff -i gene_id -u 2000 -o comparison

# Statistical analysis
hastat stat -g comparison.merged.group.csv -p phenotype.csv -o comparison.stat.csv

# Visualization
hastat plot bar --sample_hap comparison.merged.group.csv --sample_group population_groups.csv -o comparison_bar.png
```

#### 🏷️ 6. Selection sweep analysis

**Nucleotide diversity (π) analysis:**
```bash
hastat view -t pi -v test.vcf -a test.gff -i gene_id -u 2000 -g population_groups.csv -o pi_results --size 1000 --step 100
```

**Fixation index (FST) analysis:**
```bash
hastat view -t fst -v test.vcf -a test.gff -i gene_id -u 2000 -g population_groups.csv -o fst_results --size 1000 --step 100
```

#### 🏷️ 7. GWAS analysis

```bash
hastat gwas -c gwas_config.toml
```

### ⭐ Advanced Usage

**Heterozygosity filtering:**
```bash
hastat view -v test.vcf -a test.gff -i gene_id -t group --het 0.1 -o filtered_results
```

**Custom window sizes for selection analysis:**
```bash
hastat view -t pi -v test.vcf -a test.gff -i gene_id -u 1000000 -d 1000000 -g population_groups.csv -o pi_results --size 100000 --step 10000
```

**Multiple statistical test methods:**
```bash
hastat stat -g hap_groups.csv -p sample_phe.csv -m AllPairTest -o hap_stat.csv
```

**Network layout options:**
- `spring` - Spring layout (default)
- `circular` - Circular layout
- `kamada_kawai` - Kamada-Kawai layout
- `fruchterman_reingold` - Fruchterman-Reingold layout
- `shell` - Shell layout
- `spectral` - Spectral layout
- `random` - Random layout

### ⭐ Citation
If you use **hastat** in your research, please cite the following paper:

> Xiaodong Li, & Kun Lu. (2024). hastat: a python library to perform gene haplotypes analysis in natural populations. Zenodo. https://doi.org/10.5281/zenodo.11183815