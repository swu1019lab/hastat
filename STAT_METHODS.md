# Statistical repair, 2026-10-28

Based on swu1019lab/hastat commit 97736b4095d4a7b1a01fccf0371c19d2ca493a08.
This is a local repair branch, version `1.0.0-stat.2`; it does not change grouping.

On 2026-10-31, the ambiguous `AllPairTest` alias was removed at the user's
request. Both the CLI and Python API accept only `TukeyHSD` and
`WelchBonferroni`. The old name is rejected, not translated. This naming update
does not change the calculations from `1.0.0-stat.1` for the two named methods.
The CLI now writes only the requested CSV when `-o` is supplied; it does not
create a `.summary.json` sidecar. Existing historical outputs are retained.

## Models and implementation

The omnibus null is equality of all included population means; each pairwise
null is equality of two means. Independent observations and approximately normal
within-group errors are assumed. These are unadjusted one-factor comparisons.

* `TukeyHSD` (unchanged CLI default): classical equal-variance one-way ANOVA,
  `statsmodels.stats.oneway.anova_oneway(use_var='equal')`. F is the between-group
  mean square divided by the pooled within-group mean square, with k−1 and N−k
  degrees of freedom. Post-hoc Tukey–Kramer uses `pairwise_tukeyhsd` and its **raw**
  `pvalues`, `confint`, `reject` arrays; no rounded display table is parsed.
  Equal variance is required; unequal sample counts are supported by Tukey–Kramer.
* `WelchBonferroni`:
  Welch's ANOVA via explicit `use_var='unequal', welch_correction=True`.
  For w_i=n_i/s_i², W=sum(w_i), weighted mean m, define
  B=sum((1−w_i/W)²/(n_i−1)); F=sum(w_i*(mean_i−m)²)/(k−1)/
  (1+2*(k−2)*B/(k²−1)), df=(k−1,(k²−1)/(3B)).
  Two-sided Welch t tests use SciPy `ttest_ind(group2, group1, equal_var=False)`.
  For a=s1²/n1, b=s2²/n2, SE=sqrt(a+b),
  df=(a+b)²/(a²/(n1−1)+b²/(n2−1)). Differences are always **mean2−mean1**.
  For m=k*(k−1)/2 comparisons, adjusted P=min(1,m*P), and simultaneous intervals
  are difference ± t_(1−alpha/(2m),df)*SE. These Bonferroni intervals have
  family coverage at least 1−alpha under the reference model. Welch reference
  degrees of freedom are approximations, not an exact finite-sample test.

The previous implicit Welch omnibus plus pooled comparisons is replaced. There is no
automatic test-selection based on a variance pretest. Pairwise adjusted tests
are reported without an additional omnibus significance gate.

## Eligibility, edge cases and reporting

IDs must be unique on each side of a one-to-one sample join. Per-outcome finite
observations are selected **after** joining; each group must then meet the
minimum sample size (CLI default 10, website 3). Original zeros are retained;
TPM is never logged. Invalid numeric text raises an error. Missing and infinite
values are excluded and per-outcome group counts/statuses are returned.

Groups with zero sample variance are never silently removed. Welch cannot
estimate inverse-variance weights if any eligible variance is zero; return
`zero_variance`, retain group summaries, and do not fabricate P values. Classical
ANOVA retains a zero-variance group when pooled residual variance is positive;
when every variance is zero it also returns `zero_variance`. Fewer than two
eligible groups returns `insufficient_groups` and a valid empty comparison table.
Unexpected non-finite model results are `not_estimable`.

P values, effects and intervals stay in double precision; CSV uses 17 significant
digits. No small P value is rounded to zero for export. A true numerical zero
returned by the distribution implementation is retained as numerical zero, not an
invented exact significance level. The Python `result()` object records omnibus
and group status for API consumers; the CLI does not write this object to a
sidecar file.
All comparisons for one outcome form a family. Website Top 2/5/10, dragging and
renaming are display settings and do not change this family. No abc letters yet.

## References and numerical audit

* Fisher, R.A. (1925), *Statistical Methods for Research Workers*, Oliver & Boyd.
  Classical ANOVA definition: https://www.itl.nist.gov/div898/handbook/prc/section4/prc433.htm
* Tukey, J.W. (1949), Comparing Individual Means in the Analysis of Variance,
  Biometrics 5:99–114. DOI: https://doi.org/10.2307/3001913
* Kramer, C.Y. (1956), Extension of Multiple Range Tests to Group Means with
  Unequal Numbers of Replications, Biometrics 12:307–310.
  DOI: https://doi.org/10.2307/3001469
* Welch, B.L. (1951), On the Comparison of Several Mean Values: An Alternative
  Approach, Biometrika 38:330–336. DOI: https://doi.org/10.1093/biomet/38.3-4.330
* Welch, B.L. (1947), The Generalization of Student's Problem when Several
  Different Population Variances are Involved, Biometrika 34:28–35.
  DOI: https://doi.org/10.1093/biomet/34.1-2.28
* Bonferroni, C.E. (1936), Teoria statistica delle classi e calcolo delle
  probabilità, Pubblicazioni del R Istituto Superiore di Scienze Economiche e
  Commerciali di Firenze 8:3–62 (no DOI). Simultaneous interval construction:
  https://www.itl.nist.gov/div898/handbook/prc/section4/prc473.htm

`tests/test_stat.py` compares classical ANOVA with SciPy and an independently
calculated sums-of-squares formula, Tukey raw arrays with SciPy's independent
implementation, and Welch F/df/t/P/Bonferroni CIs with explicit reference formulas.
It also exercises missingness, finite-value filtering, zero variance, duplicates,
input order, empty results, tiny P values and CLI full-precision export.
`validate.py` records exact versions, source hashes, test results and checksums.
This is numerical implementation validation, not an empirical power/FDR study.
