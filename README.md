## hastat
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.11001623-blue)](https://doi.org/10.5281/zenodo.11183815)

### Description
A python library to perform gene haplotypes analysis in natural populations.

### Features
main modules:
- [x] `view` view genotype, haplotype, pi, fst, ld data of interested genes in a VCF file
- [x] `stat` perform haplotype statistic analysis with multiple phenotype using anova and multiple comparison
- [x] `plot` plot haplotype data using HapBox, HapBar, HapPie, HapNetwork, Gene, etc.
- [x] `gwas` perform gwas analysis using wrapper of gemma, plink, etc.


### Requirements
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

### Installation

```bash
git clone https://github.com/swu1019lab/hastat.git
cd hastat
# install the package using build commands (recommended)
pip install build --user
python -m build
pip install dist/hastat-0.0.5.tar.gz --user
# or
# install the package using setup.py install command
python setup.py install --user
```

### Usage
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

### Citation
If you use hastat in your research, please cite the following paper:

> Xiaodong Li, & Kun Lu. (2024). hastat: a python library to perform gene haplotypes analysis in natural populations. Zenodo. https://doi.org/10.5281/zenodo.11183815