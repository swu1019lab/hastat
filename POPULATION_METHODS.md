# VCF and population-statistics correction, 2026-10-31

This local patch follows `stat.2` and fixes the VCF reader and existing pi/Tajima's D calculations. It does not introduce a new statistical test. The package release number has not been changed. The latest complete-site window policy is validated under `qa/20261031_window_missing`, with the restored minimal CSV schema checked under `qa/20261031_minimal_output`. Earlier audits retain their historical schemas and source snapshots. Normal tool runs do not create a JSON sidecar.

Later changes: the Fst wrapper now shares the explicit window builder and has
its own validation described in [FST_METHODS.md](FST_METHODS.md). The network
repair is documented in [NETWORK_METHODS.md](NETWORK_METHODS.md). Statements
below about unchanged/out-of-scope modules describe the earlier pi/D patch.

## Coordinates and genotype calls

- All public regions, GFF gene intervals, VCF POS values and output windows use **1-based inclusive** coordinates. The numeric pysam fetch call alone uses `[start - 1, end)`. A record belongs to a region when its POS is inside it; an overlapping deletion anchored before the interval is excluded.
- Gene extensions respect strand and clip the left coordinate to 1. For a manually specified region, upstream extends toward smaller coordinates and downstream toward larger coordinates. Negative extension lengths and invalid intervals are rejected.
- The reader retains biallelic records and requires diploid called genotypes. `code1` is alternate-allele dosage, whereas `code2` retains the two allele indices. A partial call such as `0/.`, `./0`, `1/.` or `./1` becomes a whole missing genotype, consistently encoded as `-1` or `[-1, -1]`. This is a conservative software policy: the known allele of a partial call is not retained. A solitary `.` is also treated as a missing diploid call. Called haploid/polyploid genotypes and records without GT raise explicit errors.
- MAF and heterozygosity filters use completely called diploid genotypes across all VCF samples; missingness uses the number of samples as its denominator. Empty genotype results preserve the sample columns and variant index names.

## Population-statistics input and windows

Pi and Tajima's D now count chromosomes with `GenotypeArray(...).count_alleles(max_allele=1)`. Treating dosage values 0, 1 and 2 as distinct haploid alleles was incorrect. Phasing is unnecessary for these site-wise allele counts.

These two outputs use **biallelic A/C/G/T SNPs only**, excluding indels and symbolic alleles with a warning. The shared reader already excludes multiallelic records. Duplicate SNP positions are rejected, preventing split records from being counted twice. Thus these are estimates from retained biallelic SNPs, not a general estimator covering every possible variant type.

Population sample lists must be non-empty, unique, and present in the VCF. Omitting the list in the Python API uses all samples. The command-line `view` population workflow still requires a group file.

Windows begin at the requested start and advance by `step`; their inclusive end is `min(window_start + size - 1, region_end)`. Stop after the first window reaches the region end. Both size and step must be positive integers. This fixes one-base regions and prevents the last full window from acquiring an extra base through the dependency's default window builder. An overlapping-window request keeps overlap; a step larger than size deliberately leaves gaps. The existing Fst window builder is outside this patch.

`counts` means the number of input biallelic SNP records in a window after the shared reader's filters, including records monomorphic in the selected population. It is **not** the number of segregating sites used by Tajima's D, and includes incomplete SNPs subsequently excluded from D.

`--het` retains its existing filtering behavior but now warns that filtering changes the site frequency spectrum. Such results describe only the retained SNPs and should not be interpreted as an unfiltered neutrality analysis. Non-default `--maf` and `--max_missing` options are rejected for `view -t pi/tajima_d`; previously they were silently ignored. Their genotype/haplotype-output behavior is unchanged.

## Pi: estimator, denominator and missingness

For a site with allele counts c_j and n = sum(c_j), the estimated difference between two distinct sampled chromosomes is

    d = [n(n - 1) - sum_j c_j(c_j - 1)] / [n(n - 1)].

Window pi is `sum(d) / L`, where L is the inclusive window length. The implementation uses scikit-allel's `windowed_diversity` after correcting the allele-count input and supplying explicit windows. This is an estimator, with no null hypothesis, reference distribution, or p-value supplied here.

At a SNP with missing genotypes, n is the number of chromosomes from complete calls in the selected population. Estimation therefore assumes these observed calls represent that population; genotype-dependent missingness can bias it. If a retained SNP has fewer than two called chromosomes, every window containing that SNP returns NaN rather than treating unknown diversity as zero.

The VCF alone supplies no callable-base mask. The denominator assumes all requested bases are callable and that absent eligible SNP records are invariant. Consequently an empty window returns pi = 0 **under this assumption**, not as evidence of sequencing coverage. A filtered or incomplete VCF can underestimate diversity. Callable-region support remains a separate future improvement.

## Tajima's D: formula, assumptions and undefined results

Let n be the fixed number of sampled chromosomes, S the number of segregating SNPs, and K the sum of the per-site pairwise differences above. The implemented statistic is

    D = (K - S/a1) / sqrt(e1*S + e2*S*(S - 1))
    a1 = sum(i=1..n-1, 1/i); a2 = sum(i=1..n-1, 1/i^2)
    b1 = (n+1)/(3*(n-1)); b2 = 2*(n^2+n+3)/(9*n*(n-1))
    c1 = b1 - 1/a1; c2 = b2 - (n+2)/(n*a1) + a2/a1^2
    e1 = c1/a1; e2 = c2/(a1^2+a2).

For each window's retained complete sites, hastat delegates this formula to `allel.tajima_d`. The reference null is the standard neutral, constant-size, randomly mating population model with an infinite-sites mutation process. The original variance and calibration depend on the sampling and linkage model; population structure, demographic change, ascertainment and recombination affect interpretation. D is a standardized contrast, not automatically a standard-normal z score. This implementation reports D only and supplies neither calibrated p-values nor selection calls.

The selected population stays fixed across sites and windows. **Exclude each SNP with any missing diploid call among those selected samples**, then calculate the standard D from the complete sites. Missingness in samples outside the selected population does not remove a site. The number of chromosomes remains twice the selected sample count; the implementation does not estimate a variable-n correction, impute alleles, or select a different subset of individuals at each site.

This is complete-site deletion followed by the existing statistic, a common strategy discussed by Ferretti et al. (2012), not their generalized missing-data test. It avoids invalidating an entire window for one incomplete SNP. Under missingness independent of allelic state, retained sites are a less informative subset; genotype-dependent missingness or biased site ascertainment can distort the retained spectrum. A high excluded-site fraction is a limitation, not proof of adequate information simply because D is finite. No universal retained-site fraction threshold or new significance cutoff is imposed.

Windows with fewer than 3 retained segregating SNPs retain scikit-allel's minimum-sites policy. Fewer than 4 chromosomes, no complete sites, and empty windows also return NaN. The earlier policy that invalidated any incomplete window is preserved in `archive/20261031_before_window_missing_fix`, with its validation bundle. The current rule intentionally supersedes that policy.

Pi is unchanged: it can use partial coverage across individuals at a site, whereas D uses only complete sites. Thus the public window pi and D may describe different retained site sets; do not reconstruct D from the pi CSV times window length. This patch does not provide callable-base masks or missing-data-adjusted neutrality calibration.

### Tajima's D output and optional validation diagnostics

The normal CLI CSV retains exactly `chrom,start,end,counts,tajima_d,groups`. No sample-size, retained/excluded-site, status, or other auxiliary columns are added. `counts` remains the input SNP count as defined above, not S or the number of complete sites. No additional output file is produced.

The following diagnostics are available only when explicitly requesting `include_qc=True` from the Python API, as used by the developer validation scripts. They do not appear in normal CLI output:

| Column | Meaning |
|---|---|
| `n_samples` | Fixed selected diploid sample count. |
| `n_chromosomes` | Twice `n_samples`, unchanged by incomplete-site removal. |
| `sites_used` | SNP records with complete calls in the selected population, including population-monomorphic records. |
| `sites_excluded_missing` | SNP records excluded because at least one selected sample is missing; `counts = sites_used + sites_excluded_missing`. |
| `segregating_sites` | S in the D formula, computed only among `sites_used`. |

The Python methods `get_tajima_d_data` and `get_multi_loci_tajima_d_data` retain their five-column default result. Set `include_qc=True` for columns defined by `TAJIMA_QC_COLUMNS` only when those diagnostics are specifically needed. The CLI uses the five-column default plus its existing groups column. Undefined statistics remain NaN. A finite D is not a biological quality certification or a significant test result. Execution and validation status belong in developer QA records, not normal tool result files.

The production D values come from `GeneVariant.get_tajima_d_data` calling `allel.tajima_d` on the selected complete-site allele counts. Independent formulas in the validation scripts are comparison references only; they do not supply or replace tool output values.

## References and implementation mapping

1. Nei M, Li WH (1979). *Mathematical model for studying genetic variation in terms of restriction endonucleases*. PNAS 76:5269–5273. DOI: [10.1073/pnas.76.10.5269](https://doi.org/10.1073/pnas.76.10.5269). Basis for nucleotide diversity as average pairwise differences; the finite-sample pair-count expression above is evaluated over distinct chromosomes.
2. Tajima F (1989). *Statistical method for testing the neutral mutation hypothesis by DNA polymorphism*. Genetics 123:585–595. DOI: [10.1093/genetics/123.3.585](https://doi.org/10.1093/genetics/123.3.585). Basis for the standardized contrast and variance constants above.
3. [scikit-allel diversity implementation](https://scikit-allel.readthedocs.io/en/stable/_modules/allel/stats/diversity.html): `mean_pairwise_difference`, `windowed_diversity`, `tajima_d`; validation pins scikit-allel 1.3.13 and records the installed implementation hashes.
4. [pysam VariantFile.fetch](https://pysam.readthedocs.io/en/latest/api.html#pysam.VariantFile.fetch): numeric fetch coordinates are zero-based half-open. [scikit-allel window implementation](https://scikit-allel.readthedocs.io/en/stable/_modules/allel/stats/window.html) motivates supplying explicit inclusive windows.
5. Ferretti L, Raineri E, Ramos-Onsins S (2012). *Neutrality Tests for Sequences with Missing Data*. Genetics 191:1397–1401. DOI: [10.1534/genetics.112.139949](https://doi.org/10.1534/genetics.112.139949). Discusses removing incomplete positions/individuals and its information loss, and derives more general missing-data statistics. This patch implements only complete-site filtering before standard Tajima D, not the generalized formulas from that article.

## Validation and scope

`tests/test_population.py` builds real bgzip/tabix VCF fixtures. Pi is checked by enumerating chromosome pairs for all 27 three-sample diploid genotype combinations; Tajima's D is checked against an independent formula across 20 synthetic matrices and the library's correctly counted input. Tests also cover boundaries, strand, partial calls, empty/short windows, sample selection, CSV output and a shared-reader Fst regression. Existing stat tests run in the same validation session.

The latest policy also tests partial-call removal, selected-population specificity, fixed sample membership, retained/excluded counts, overlapping windows and 15 synthetic incomplete genotype matrices against independent complete-site formulas. `validate_window_missing.py --fixtures-dir <hastat-window-comparison-20261031/inputs>` runs the suite and repeats the same frozen 780 windows at 1/5/10 kb, checking that pi and complete-data D are unchanged. Execute it only inside an authorized PBS core/A6000 allocation. The validation-only output directory records strict-policy/new-policy comparisons, test status, input/source hashes and environment versions. Numerical agreement verifies implementation; it is not an empirical error-rate or biological validation study. No real sample data are used. Network algorithms, haplotype grouping semantics and the Fst estimator are unchanged. Do not rerun the earlier `validate_population.py` into its frozen historical output directory to validate this new policy.
