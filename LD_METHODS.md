# LD matrix repair (2026-10-31)

## Estimator and interpretation

Production off-diagonal values continue to use scikit-allel 1.3.13
`allel.rogers_huff_r`, not a newly introduced test. Biallelic diploid genotypes
are coded as alternate-allele dosage 0/1/2, with -1 for missing. The existing VCF
reader treats a partially called genotype as wholly missing. Biallelic indels
remain included; multiallelic records remain excluded by the reader.

For each pair, retain only samples called at both sites. On those samples,
r = sum((x-xbar)*(y-ybar)) / sqrt(sum((x-xbar)^2)*sum((y-ybar)^2)); return r².
The library uses the algebraically equivalent moment-sum formula. There is no
imputation, no global complete-sample filter, and no replacement of NaN by zero.
Different pairs can therefore use different samples. With two shared variable
observations r² is mathematically 1; this is not evidence of a reliable population
estimate. No P values, null distribution, significance test, confidence intervals
or correction for population structure/relatedness are introduced here.

This is an unphased genotype-based estimate, not direct counting of phased
haplotypes. Biological interpretation as gametic LD requires the estimator's
population-model conditions; do not equate genotype correlation across structured
or related samples with physical linkage. Cross-chromosome matrices are permitted,
but correlations between chromosomes do not establish physical linkage.

Reference: Rogers AR & Huff C (2009), Linkage Disequilibrium Between Loci With
Unknown Phase. Genetics 182:839-844. DOI: https://doi.org/10.1534/genetics.108.093153.
The library docstring labels the method 2008; the published article is 2009.
API and pinned implementation:
https://scikit-allel.readthedocs.io/en/stable/stats/ld.html
https://github.com/cggh/scikit-allel/blob/v1.3.13/allel/stats/ld.py
https://github.com/cggh/scikit-allel/blob/v1.3.13/allel/opt/stats.pyx

## Boundary cases and output

- Both single/multi-region APIs return an n_variants x n_variants NumPy array.
  An empty region is (0,0), not (1,1) or a one-dimensional empty vector.
- Diagonal entries are 1 only if at least two called diploid samples have
  different dosages at that site. All-reference, all-alternate, all-heterozygous,
  all-missing or only-one-called-sample sites have NaN diagonals.
- For an off-diagonal pair, insufficient shared observations or constant dosage
  at either site within the shared subset leaves NaN, even if both diagonals are 1.
- Nonfinite library results remain NaN. Finite r² values are clipped to [0,1] to
  contain float32 roundoff (the documented library example exceeds 1 slightly).
  For ordinary finite off-diagonal values <=1, the existing float32 calculation
  and output precision are retained exactly. Returning float64 would not recover
  precision lost inside the library, so no such precision claim is made.
- Single-region order is VCF record order. Multi-region requests use
  normalize_loci: chromosome names sorted lexicographically, intervals ordered
  and overlaps/adjacency merged; VCF tie order is retained within each interval.
  Overlapping requests do not duplicate records. Exact duplicate full
  (chrom,pos,ref,alt) records in the VCF raise an error.
- Calling the multi-region API with no requested intervals is an input error;
  valid requested intervals containing no variants return (0,0).

LD remains a Python API feature. There is no `view --type ld` CLI option;
this patch does not add an output file type, status column or JSON sidecar.

```python
import pandas as pd
from hastat.utils.haplotype import normalize_loci
from hastat.viz import ldheatmap

# gv is a GeneVariant backed by an indexed pysam.VariantFile.
loci = normalize_loci([('chr1', 1, 1000), ('chr1', 800, 2000)])
genotypes = gv.get_multi_loci_geno_data(loci)
r2 = gv.get_multi_loci_ld_data(loci)
frame = pd.DataFrame(r2, index=genotypes.index, columns=genotypes.index)
ax = ldheatmap.plot(frame, plot_value=True, plot_snp=True)
ax.figure.savefig('ld.pdf', bbox_inches='tight')
```

## Heatmap behavior and validation

The heatmap requires a nonempty square symmetric numeric DataFrame with matching
unique row/column labels; values are in [0,1] or NaN. Undefined entries are gray,
with a legend, and do not receive numeric annotations. A single variant can show
its diagonal; asking to hide its only cell raises an explicit error. A physical
position track requires one chromosome; a single or repeated position is centered
instead of dividing by zero. Colors for defined r² retain the existing colormap.

`tests/test_ld.py` enumerates all 256 length-four vectors over -1/0/1/2,
checking the resulting 65,536 ordered matrix entries against an independent
double-precision centered-dot-product reference, with explicit NaN masks.
Additional indexed VCF tests cover coordinates, partial calls, zero/one/two
variants, overlapping requests, duplicate records, indels and chromosome order.
Real library call interception verifies the production estimator remains in use.
Rendering tests cover gray undefined cells, all-NaN input and position tracks.
`validate_ld.py` compares the archived code with this patch on fixed synthetic
inputs and runs all prior regressions in a PBS core/A6000 allocation.
QA records and example exports are developer-only, not normal user outputs.
