# Plot data and statistical annotations (2026-10-31)

## Inputs and outputs

Box/bar/pie CSVs use the first column as the sample ID. Haplotype input has
exactly two columns; phenotype/group metadata has at least two. CSV fields
are read as strings, so 001, 1 and NA remain distinct identifiers. Surrounding
whitespace is stripped from sample and category labels. Empty/missing labels,
duplicate IDs (including duplicates created by trimming), duplicate headers
and reserved metadata column names sample/haplotypes are rejected. Validation
precedes haplotype display filtering. DataFrame inputs are copied, not mutated.

The join is one-to-one on shared samples; valid subsets are allowed and numbers
of unmatched rows are logged. No common samples or no finite displayed phenotype
values is an error. Only phenotype values are converted to numeric; conventional
missing tokens and nonfinite values are excluded per trait, zeros retained,
and other nonnumeric text raises an error. Bar/pie values are counts or the
existing within-population proportions of the retained plotted samples.

Normal output remains only the requested figure. There are no additional status
columns, tables or JSON sidecars. Validation evidence is confined to developer QA.
Plot failures propagate to a nonzero CLI exit instead of being logged as success.

## Comparison families and plotting

The default remains t-test. Tests are independent-sample, two-sided tests.
The family is every pair among k groups with at least --min_hap_size finite
matched observations for the selected phenotype, m=k(k-1)/2. The plot default
minimum is 2; stat defaults to 10. Comparisons and haplotypes select displayed
brackets/boxes only, so hiding a group does not reduce m or alter its P values.
Tests across different traits or separate invocations are not jointly corrected.
Changing the supplied input population can change the inference.

TukeyHSD and WelchBonferroni call the existing HapAnovaTest class with the same
observations and minimum. See STAT_METHODS.md for original references, formulas,
assumptions and numerical validation. They use its full precision p-adj array,
without rounded tables or a second correction. The existing welch name delegates
to WelchBonferroni for plotting, including its zero-variance family guard.
Equality to stat requires identical observations, groups, method and minimum.
Pooled pairwise t tests are not Tukey tests when more than two groups are used.

The legacy t-test and mannwhitneyu methods now use min(1,m*p_raw). Unestimable
pairs do not shrink m. Two-sided pooled t uses
s_p^2=((n1-1)s1^2+(n2-1)s2^2)/(n1+n2-2),
t=(mean1-mean2)/sqrt(s_p^2*(1/n1+1/n2)), df=n1+n2-2.
The null is equal means; observations must be independent, with normal errors
and equal population variances for the t reference distribution. The implementation
is scipy.stats.ttest_ind(equal_var=True, alternative='two-sided'). Both group
variances zero, too few samples or nonfinite results yield no annotation.

Mann-Whitney tests equality of the two underlying distributions; interpreting it
as only a median difference requires additional shape assumptions. Observations
are independent and ordinal/comparable. U=sum[ I(x>y)+0.5*I(x=y) ]. The existing
SciPy auto rule is made explicit: exact U null distribution when either n<=8
and there are no ties; otherwise a tie-adjusted normal approximation with 0.5
continuity correction. Its variance is n1*n2/12 *
[N+1-sum(t_j^3-t_j)/(N*(N-1))]. Small tied samples retain this approximation;
it is not an exact permutation test and no new calibration claim is made.

Bonferroni controls family-wise error by the union bound when component P values
are valid; it is not an FDR adjustment. Numerical tests compare t with an explicit
formula, Mann-Whitney with enumerated exact U and a tie-corrected normal formula,
Tukey with SciPy's separate implementation and both named methods with stat.
This validates implementation, not empirical error control in structured populations.

Stars use adjusted P: *** <0.001, ** <0.01, * <0.05, otherwise ns. NaN, infinity,
insufficient samples and non-estimable families are never ns: skip the bracket
and warn in the log. Small groups remain visible as observations. Unknown pairs,
self comparisons and repeated unordered pairs are rejected. Brackets are placed
above the data even for negative outcomes. Jitter is deterministic and does not
reset the caller's global NumPy random state.

## References

- Student (1908), The Probable Error of a Mean. Biometrika 6:1-25.
  https://doi.org/10.1093/biomet/6.1.1 (Student t reference distribution).
- Mann & Whitney (1947), On a Test of Whether one of Two Random Variables is
  Stochastically Larger than the Other. Annals of Mathematical Statistics 18:50-60.
  https://doi.org/10.1214/aoms/1177730491
- Bonferroni (1936), Teoria statistica delle classi e calcolo delle probabilita,
  Pubblicazioni del R Istituto Superiore di Scienze Economiche e Commerciali di
  Firenze 8:3-62 (no DOI); see STAT_METHODS.md for Tukey/Kramer/Welch references.
- Pinned API behavior:
  https://docs.scipy.org/doc/scipy-1.15.2/reference/generated/scipy.stats.ttest_ind.html
  https://docs.scipy.org/doc/scipy-1.15.2/reference/generated/scipy.stats.mannwhitneyu.html

Reproduce the synthetic audit with validate_plot.py in a core/A6000 PBS allocation.
QA figures show native hastat output for regression inspection, not manuscript
panels. Existing tool colors are retained; no OmicsX main figure is changed.

## Vector text export

All plot savers use a shared, scoped Matplotlib rc_context with
`pdf.fonttype=42` and `svg.fonttype=none`. PDF embeds TrueType subsets
(Type 0/CIDFontType2 in the PDF structure); SVG retains text elements.
Caller settings are restored even if saving fails. PNG output is unchanged.
No new CLI argument, table column or result-side summary file is added.

These settings preserve text for downstream editing, but Illustrator still
needs the fonts installed locally or substituted. Mathematical and rotated
labels can be split across text objects. Validation checks embedded fonts,
text extraction and SVG text nodes; it does not assert a completed editing
round-trip inside Illustrator.
