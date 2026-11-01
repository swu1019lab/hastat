# Weir-Cockerham Fst wrapper correction, 2026-10-31

This patch changes input validation, window boundaries and treatment of
unestimable records. It retains Weir and Cockerham's theta estimator and calls
`allel.windowed_weir_cockerham_fst` from scikit-allel (validation: 1.3.13).
It adds no statistical test, P value, significance threshold or result column.

## Estimator and implementation

For each biallelic locus and allele, let r=2, n_i be the number of completely
called diploid individuals in population i, p_i its allele frequency,
n_bar=sum(n_i)/r, p_bar=sum(n_i*p_i)/sum(n_i), and h_bar the fraction of called
individuals heterozygous for the allele. Then:

```
n_C = [sum(n_i) - sum(n_i^2)/sum(n_i)] / (r-1)
s2  = sum[n_i*(p_i-p_bar)^2] / [(r-1)*n_bar]
a   = (n_bar/n_C) * {s2 - [p_bar*(1-p_bar) - (r-1)*s2/r - h_bar/4]/(n_bar-1)}
b   = n_bar/(n_bar-1) * {p_bar*(1-p_bar) - (r-1)*s2/r - (2*n_bar-1)*h_bar/(4*n_bar)}
c   = h_bar/2
theta_window = sum(a) / [sum(a) + sum(b) + sum(c)]
```

Sums extend over retained loci and alleles. Window Fst is a ratio of component
sums, not the arithmetic mean of per-site ratios. Negative estimates are valid
sampling outcomes and remain unmodified; non-finite results become NaN.
This is a diploid population-structure estimator, not a direct selection test.
No null reference distribution or P value is asserted. Sampling, relatedness,
population definitions, genotype error and non-random missingness affect its
interpretation. Linkage matters for uncertainty, which this output does not
estimate.

The production calculation remains in scikit-allel. The scalar formulas in
`tests/test_fst.py` are independent QA references and never supply normal output.

## Missing data and sample validation

Sample lists must be non-empty, unique, present in the VCF and disjoint.
The CLI group file must contain two columns and at least two populations; each
sample appears exactly once. Literal names such as `001` and `NA` are preserved.
Malformed assignments cause an error rather than a silently smaller population.

The shared reader treats a partial call (e.g. 0/.) as wholly missing. Each site
uses the remaining complete diploid calls **within the two selected populations**.
This is not the complete-site policy used for Tajima's D.

For r=2, both populations must have at least one called individual, and more
than two called individuals must remain in total, ensuring n_bar>1 and n_C>0.
Sites failing these algebraic prerequisites are excluded from all a/b/c sums
before calling the library. One individual versus two can be estimable; one
versus one cannot. No arbitrary universal two-per-population cutoff is imposed.

This exclusion matters: the library sums non-NaN a, b and c separately. At an
unestimable site, a and b can be NaN while c is finite. Passing such sites
unfiltered can change the denominator of a mixed window or report a misleading
zero for an unestimable window. The wrapper now filters these records together.
An empty window, a window with no estimable records or with zero total
denominator returns NaN. It does not become zero differentiation by default.

Missingness need not invalidate an entire window, but the estimate describes
the observed and estimable records; genotype-dependent missingness can bias it.
The number of excluded records is logged, not added to the normal CSV.

## Records, windows and output

The existing biallelic-record scope is unchanged: SNPs and biallelic indels can
contribute. The shared reader excludes multiallelic records and requires diploid
GT. No additional SNP-only conversion or variant deduplication is introduced.
`--het` retains its existing filter across VCF samples, so estimates then
describe that filtered marker set. Non-default `--maf` and `--max_missing` are
explicitly rejected for Fst, as they already are for pi/D; these parameters
previously had no effect for Fst. This patch does not invent a new filtering
policy by activating them silently.

All coordinates are 1-based inclusive. The shared explicit window builder uses
`end=min(start+size-1, region_end)` and stops when a window reaches region_end.
It handles a one-base region, a one-base final window, overlap and deliberately
gapped steps. Size and step must be positive integers. Separate loci are
processed separately and stacked with their actual chromosome/coordinates.

Normal CLI output is exactly `chrom,start,end,counts,fst,groups`. `counts`
remains the number of input biallelic records in the window after shared-reader
filters, including unestimable records. It is not an effective-site count.
No status column, summary sidecar, or additional normal output is generated.
NumPy error handling is scoped with `errstate`, including error paths.

## References and validation

- Weir BS, Cockerham CC (1984). Estimating F-statistics for the analysis of
  population structure. Evolution 38:1358–1370.
  [DOI: 10.1111/j.1558-5646.1984.tb05657.x](https://doi.org/10.1111/j.1558-5646.1984.tb05657.x).
- [Official scikit-allel F-statistics documentation](https://scikit-allel.readthedocs.io/en/stable/stats/fst.html).
- [Official implementation and component mapping](https://scikit-allel.readthedocs.io/en/stable/_modules/allel/stats/fst.html).
  Validation also records hashes of the installed 1.3.13 fst/window source.

`validate_fst.py` executes `tests/test_fst.py` plus existing stat/population/network
tests inside PBS core/A6000. It compares fixed synthetic inputs at 1/5/10 kb
against the archived wrapper and independent scalar formulas. See
`qa/20261031_fst/PLAN.md` for the frozen validation conditions; output and source
hashes are developer QA, not routine tool outputs. This is numerical validation,
not a biological benchmark or calibration of an empirical error rate.
