## hastat

### Description
An easy library to perform gene haplotypes analysis in Python.

### Features
- [x] Generate haplotypes from genotypes
- [x] Calculate haplotype frequencies
- [x] Calculate linkage disequilibrium and Fst
- [x] Calculate haplotype association with a phenotype using ANOVA

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
pip install dist/hastat-0.0.3.tar.gz --user
# or
# install the package using setup.py install command
python setup.py install --user
```

### Usage

```bash
hastat --ann <gff_file> --gen <vcf_file> --phe <phenotype_file> --out-dir <output_dir> genes.list
```

