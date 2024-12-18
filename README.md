![hastat](tests/hastat.png)

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.11001623-blue)](https://doi.org/10.5281/zenodo.11183815)

## Table of Contents
- [hastat](#hastat)
  - [Description](#-description)
  - [Features](#-features)
  - [Requirements](#-requirements)
  - [Installation](#-installation)
  - [Usage](#-usage)
  - [Example](#-example)
    - [1. view the genotype or haplotypes data of a gene](#-1-view-the-genotype-or-haplotypes-data-of-a-gene)
    - [2. perform haplotype statistic analysis](#-2-perform-haplotype-statistic-analysis)
    - [3. plot the haplotype data of a gene](#-3-plot-the-haplotype-data-of-a-gene)
      - [bar plot](#-1%EF%B8%8F%E2%83%A3-bar-plot)
      - [pie plot](#-2%EF%B8%8F%E2%83%A3-pie-plot)
      - [box plot](#-3%EF%B8%8F%E2%83%A3-box-plot)
      - [gene plot](#-4%EF%B8%8F%E2%83%A3-gene-plot)
  - [Citation](#-citation)

### ⭐ Description
A python library to perform gene haplotypes analysis in natural populations.

### ⭐ Features
main modules:
- [x] `view` view genotype, haplotype, pi, fst, ld data of interested genes in a VCF file
- [x] `stat` perform haplotype statistic analysis with multiple phenotype using anova and multiple comparison
- [x] `plot` plot haplotype data using HapBox, HapBar, HapPie, HapNetwork, HapGene, etc.
- [x] `gwas` perform gwas analysis using wrapper of gemma, plink, etc.


### ⭐ Requirements
Python 3.9 or higher and the following packages are required:
- pandas
- numpy
- scipy
- statsmodels
- scikit-allel
- gffutils
- pysam
- matplotlib
- tomli
- prettytable

### ⭐ Installation

```bash
git clone https://github.com/swu1019lab/hastat.git
cd hastat
# install the package using build commands (recommended)
pip install build --user
python -m build
pip install dist/hastat-0.0.6.tar.gz --user
# or
# install the package using setup.py install command
python setup.py install --user
```

### ⭐ Usage
```bash
usage: hastat [-h] [--version] [--log LOG] {view,stat,plot,gwas} ...

A package for gene haplotype analysis

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --log LOG             The log file name (default: stdout)

subcommands:
  valid subcommands

  {view,stat,plot,gwas}
                        additional help
    view                View the haplotypes data of a gene or target region
    stat                Perform statistical analysis on the haplotypes of a gene
    plot                Plot the haplotypes data of a gene
    gwas                Perform GWAS analysis using GEMMA/EMAX wrapper
```

### ⭐ Example
#### 🏷️ 1. view the genotype or haplotypes data of a gene
```bash
hastat view -v test.vcf -i gene_id -t hap_group -o hap_groups.csv
```
the `hap_groups.csv` file will contain the haplotype data of the gene with the following format:
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
other options can be found by running `hastat view -h`:

```bash
usage: hastat view [-h] -v VCF [-a GFF] (-r REGION | -i GENE_ID) [-t {genotype,hap_table,hap_group,hap_freq}] [-g GROUP] [-o OUT]

optional arguments:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     The VCF file containing the genotype data
  -a GFF, --gff GFF     The GFF file containing the gene annotation
  -r REGION, --region REGION
                        The region of the gene to be analyzed, format: chr:start-end
  -i GENE_ID, --gene_id GENE_ID
                        The gene ID for the target region
  -t {genotype,hap_table,hap_group,hap_freq}, --type {genotype,hap_table,hap_group,hap_freq}
                        The data type to be analyzed (default: hap_group)
  -g GROUP, --group GROUP
                        A csv file containing the custom groups of samples if the data type is hap_group
  -o OUT, --out OUT     The output csv file name (default: stdout)
```
#### 🏷️ 2. perform haplotype statistic analysis
```bash
hastat stat -g hap_groups.csv -p sample_phe.csv -o hap_stat.csv
```
the `hap_groups.csv` file will contain the haplotype groups of the gene with the following format:
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
the `sample_phe.csv` file will contain the phenotype data of the gene with the following format:
```csv
samples,trait
sample1,5
sample2,3
sample3,7
sample4,2
sample5,9
sample6,4
sample7,6
sample8,8
```

the `hap_stat.csv` file will contain the haplotype statistic data of the gene with the following format:
```csv
group1,group2,meandiff,p-adj,lower,upper,reject,pheno,annotate,anova,count1,count2,mean1,mean2
Hap1,Hap2,2.0,0.7269,-3.4809,7.4809,False,trait,gene,0.6503121273263954,5.0,5.0,4.2,6.2
Hap1,Hap3,0.8,0.9747,-4.6809,6.2809,False,trait,gene,0.6503121273263954,5.0,5.0,4.2,5.0
Hap1,Hap4,2.4,0.6041,-3.0809,7.8809,False,trait,gene,0.6503121273263954,5.0,5.0,4.2,6.6
Hap2,Hap3,-1.2,0.9221,-6.6809,4.2809,False,trait,gene,0.6503121273263954,5.0,5.0,6.2,5.0
Hap2,Hap4,0.4,0.9966,-5.0809,5.8809,False,trait,gene,0.6503121273263954,5.0,5.0,6.2,6.6
Hap3,Hap4,1.6,0.837,-3.8809,7.0809,False,trait,gene,0.6503121273263954,5.0,5.0,5.0,6.6
```
other options can be found by running `hastat stat -h`:
```bash
usage: hastat stat [-h] -g GROUP -p PHENO [-s MIN_HAP_SIZE] [-a ANNOTATE] [-m {TukeyHSD,AllPairTest}] [-o OUT]

optional arguments:
  -h, --help            show this help message and exit
  -g GROUP, --group GROUP
                        A csv file containing the haplotype groups of samples
  -p PHENO, --pheno PHENO
                        A csv file containing the phenotype data of samples
  -s MIN_HAP_SIZE, --min_hap_size MIN_HAP_SIZE
                        The minimum sample size for each haplotype (default: 10)
  -a ANNOTATE, --annotate ANNOTATE
                        The annotation of haplotypes (default: gene)
  -m {TukeyHSD,AllPairTest}, --method {TukeyHSD,AllPairTest}
                        The method for multiple comparisons (default: TukeyHSD)
  -o OUT, --out OUT     The output csv file name (default: stdout)
```
#### 🏷️ 3. plot the haplotype data of a gene
##### 1️⃣ bar plot
`bar.toml` file:
```toml
[data]
sample_hap = 'hap_groups.csv' # a csv file with the haplotype frequencies
sample_group = 'sample_groups.csv' # a csv file with the sample group information

[plot]
group_index = 1 # 0-based index of the column in the sample group file that contains the group names
save_fig = 'hastat_bar.png' # can be 'hastat_bar.png' or 'hastat_bar.pdf'
width = 10 # set the width of the figure in inches
height = 5 # set the height of the figure in inches
color = ['#498DCB', '#F9BEBF', '#747474', '#EE3424'] # color for the bars
x_label = 'Haplotype'
y_label = 'Haplotype frequency'
calc_percentage = true # calculate the percentage of each haplotype
```

```bash
hastat plot -t bar -c bar.toml
```
<div align="center">
  <img alt="hastat_bar.png" height="200px" src="tests/hastat_bar.png"/>
</div>

##### ️2️⃣ pie plot
`pie.toml` file:
```toml
[data]
sample_hap = 'hap_groups.csv' # a csv file with the haplotype frequencies
sample_group = 'sample_groups.csv' # a csv file with the sample group information

[plot]
group_index = 1 # 0-based index of the column in the sample group file that contains the group names
save_fig = 'hastat_pie.png' # can be 'hastat_bar.png' or 'hastat_bar.pdf'
width = 5 # set the width of the figure in inches
height = 2 # set the height of the figure in inches
calc_percentage = false # calculate the percentage of each haplotype
```

```bash
hastat plot -t pie -c pie.toml
```
<div align="center">
  <img alt="hastat_pie.png" height="100px" src="tests/hastat_pie.png"/>
</div>

##### ️3️⃣ box plot
`box.toml` file:
```toml
[data]
sample_hap = 'hap_groups.csv'
sample_phe = 'sample_phe.csv'

[plot]
group_index = 1 # 0-based index of the column in the sample phenotype file that contains the trait names
comparisons = [['Hap1', 'Hap2'], ['Hap1', 'Hap3'], ['Hap3', 'Hap4']]
sig_symbol = ['*', '**', '***']
step_size = [0.1, 1.5, 0.1]
save_fig = 'hastat_box.png'
width = 5
height = 5
```

```bash
hastat plot -t box -c box.toml
```
<div align="center">
  <img alt="hastat_box.png" height="200px" src="tests/hastat_box.png"/>
</div>

##### ️4️⃣ gene plot
`gene.toml` file:
```toml
[data]
ann_file = "gene_ann.csv"
hap_file = "gene_hap.csv"

[plot]
show_hap = ['Hap1', 'Hap2', 'Hap3', 'Hap4']
name = 'gene1'
chrom = 'Chr01'
start = 9483087
end = 9486387
strand = '+'
save_fig = 'hastat_gene.png'
width = 8
height = 8
```

`gene_ann.csv` file:
```csv
chrom,start,end,strand,name,type
Chr01,9483087,9484131,+,gene1,exon
Chr01,9483087,9484131,+,gene1,CDS
Chr01,9484131,9484219,+,gene1,intron
Chr01,9484219,9485730,+,gene1,exon
Chr01,9484219,9485730,+,gene1,CDS
Chr01,9485730,9485801,+,gene1,intron
Chr01,9485801,9486387,+,gene1,exon
Chr01,9485801,9486387,+,gene1,CDS
```

`gene_hap.csv` file:
```csv
chrom,Chr01,Chr01,Chr01,Chr01,haplotypes
pos,9481937,9482027,9483706,9484042,
ref,T,A,A,C,
alt,C,G,G,A,
samples,,,,,
sample1,0,0,0,0,Hap1
sample2,2,2,0,0,Hap2
sample3,0,2,2,0,Hap3
sample4,0,1,1,0,Hap4
sample5,0,0,0,0,Hap1
sample6,2,2,0,0,Hap2
sample7,0,2,2,0,Hap3
sample8,0,1,1,0,Hap4
...
```

```bash
hastat plot -t gene -c gene.toml
```
<div align="center">
  <img alt="hastat_gene.png" height="200px" src="tests/hastat_gene.png"/>
</div>

other options can be found by running `hastat plot -h`:
```bash
usage: hastat plot [-h] -t {bar,pie,box,gene,network} -c CONFIG

optional arguments:
  -h, --help            show this help message and exit
  -t {bar,pie,box,gene,network}, --type {bar,pie,box,gene,network}
                        The plot type
  -c CONFIG, --config CONFIG
                        The configuration file for plotting
```

#### 🏷️ 4. perform GWAS analysis


### ⭐ Citation
If you use **hastat** in your research, please cite the following paper:

> Xiaodong Li, & Kun Lu. (2024). hastat: a python library to perform gene haplotypes analysis in natural populations. Zenodo. https://doi.org/10.5281/zenodo.11183815