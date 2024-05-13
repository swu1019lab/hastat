## hastat

### Description
A python library to perform gene haplotypes analysis in natural populations. It can generate haplotypes from genotypes, calculate haplotype frequencies, linkage disequilibrium, Fst, and haplotype association with a phenotype using ANOVA. It also provides visualisation tools to display the analysis results, such as HapBox, HapBar, Gene, etc. The library is designed to be user-friendly and easy to use, and it can be easily integrated into existing pipelines for gene haplotypes analysis in natural populations.

### Features
- [x] Generate haplotypes from genotypes (VCF file)
- [x] Calculate haplotype frequencies
- [x] Calculate linkage disequilibrium and Fst
- [x] Calculate haplotype association with a phenotype using ANOVA
- [x] Visualise analysis results using HapBox, HapBar, Gene, etc.

### Requirements
- Python 3.9 or higher
- pandas
- numpy
- scipy
- statsmodels
- scikit-allel
- gffutils
- pysam

### Installation

```bash
git clone https://github.com/swu1019lab/hastat.git
cd hastat
# install the package using build and pip commands
pip install build --user
python -m build
pip install dist/hastat-0.0.4.tar.gz --user
# or
# install the package using setup.py install command
python setup.py install --user
```

### Usage

```bash
hastat --ann <gff_file> --gen <vcf_file> --phe <phenotype_file> --out-dir <output_dir> genes.list
```

