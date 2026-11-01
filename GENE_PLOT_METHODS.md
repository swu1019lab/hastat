# Gene plot inputs, clustering and window tracks

## Sample counts on haplotype rows

Set `hap_file = "gene"` in the gene TOML. The plot reads `gene.table.csv` and
discovers the matching `gene.group.csv` automatically. An explicit
`hap_file = "gene.table.csv"` also derives this same prefix. Counts from the
`view -t group` sample-to-haplotype CSV appear as `Hap1 (n=25)`;
`hap_count_show = false` hides them. If the group file is absent, labels retain
their original names, including when `hap_count_show = true`. A missing table
or a present but unreadable/malformed group file is an error. Present group files
are validated even if counts are hidden. No directory scan or cross-prefix
fallback is used. The former `sample_hap_file` setting is no longer accepted.
`hap_label_show = false` hides the entire label, including the count.

Count one row per unique sample after validating literal sample/haplotype IDs;
duplicate sample IDs are rejected, including duplicates outside a displayed
subset. These are sample counts for the supplied assignments, not chromosome
copy counts, SNP counts, or phenotype-test sample sizes. Use table/group files
generated with the same gene/region, VCF sample set and filtering parameters.

Selecting, filtering, grouping or clustering displayed haplotypes does not merge
their assignments or recompute their counts. Assignments to hidden haplotypes are
ignored for display. Every displayed haplotype must have at least one assignment;
an incomplete or mismatched map raises an error instead of inventing zero counts.
Longer labels are measured before export; group bars and clustering trees are
moved left as needed without shrinking the matrix. No result CSV columns or
automatic sidecars are added. Assignments to a haplotype absent from the entire
input table are rejected, even if it would be hidden. Table/group IDs alone cannot
prove common provenance: users must generate the files from the same inputs.

## Category colors

`hap_color_file` replaces the former `hap_group_file` setting. It is an optional,
headerless, exactly two-column haplotype-to-category CSV, e.g. `Hap1,A` / `Hap2,B`.
Colors are assigned to categories automatically. Empty fields and duplicate
haplotype IDs (even repeated identical assignments) are rejected. Existing
category selection and ordering remain unchanged. Old parameter names raise
actionable migration errors rather than being ignored.

## Haplotype display clustering

Enable `hap_cluster_show = true`. The default `hap_cluster_linkage` is `average`
(UPGMA); `single`, `complete` and `weighted` are also supported. These operate
directly on dissimilarities. `ward`, `centroid` and `median` are rejected because
the supplied Hamming distances are not guaranteed Euclidean. Previously Ward
was the default despite using precomputed Hamming distances.

For the displayed profiles, remove every locus missing in any displayed
haplotype, and calculate d(i,j) = number of unequal genotype states / number of
common complete loci. States 0, 1 and 2 are discrete categories; -1 and explicit
missing tokens are missing, not a fourth biological state. Average linkage joins
clusters using the mean of all cross-cluster pairwise distances. Every distinct
profile has equal weight; the number of assigned samples affects labels only.
The table reader rejects invalid states, duplicate variant records/headers and
invalid coordinates before plotting.

IDs are sorted before distance calculation to make ties reproducible under input
column permutations for a fixed SciPy version; equally valid tied trees can still
differ across implementations. A two-haplotype input produces one branch; a
singleton has no tree. Identical displayed profiles remain separate, with zero
distance. With two or more profiles, no common complete loci causes an error.
Invalid methods/errors are not swallowed as a successful plot without a tree.
Labels, counts and category bars follow the final leaf order, with tree spacing
measured for both name-only and count labels.

This is descriptive genotype-profile clustering, not a phylogeny, ancestry
inference, bootstrap support or a newly introduced significance test. No null
hypothesis, reference distribution or p-value is claimed.

Implementation reference: [SciPy 1.15.2 linkage documentation](https://docs.scipy.org/doc/scipy-1.15.2/reference/generated/scipy.cluster.hierarchy.linkage.html),
including the average-linkage formula and the Euclidean-distance restriction.
Algorithm reference cited by SciPy: D. Mullner, *Modern hierarchical,
agglomerative clustering algorithms* (2011), [arXiv:1109.2378](https://arxiv.org/abs/1109.2378).
Numerical validation uses hand-calculated Hamming distances 1/4, 1, 3/4 and
average-linkage heights 1/4, 7/8, plus tie permutations, missing loci and edge cases.

## Shared physical coordinates

Gene structures, annotation targets, physical SNP connection endpoints and
Pi/Fst tracks use the same selected interval, including strand-aware upstream
and downstream extensions. Physical coordinates always increase left to right;
only the strand arrow changes direction. For inclusive interval [start, end],
base centers are integers and the outer edges are start - 0.5 and end + 0.5.
The drawing position is x0 + width * (position - start + 0.5) / (end-start+1).
Features cover [feature_start-0.5, feature_end+0.5], preserving inclusive length
and making single-base genes/features valid. Gene names center over the gene
body. Annotation targets in the selected flanks are also retained.
Strand direction heads are placed inside uncovered gaps between merged feature
blocks. Gaps too narrow for a readable head are skipped. When no readable gap
exists, a short horizontal direction marker appears above the structure.
Direction markers do not change feature coordinates or cover feature blocks.

The physical frame starts at bbox.x + xy[0] * bbox.width and spans
width * bbox.width. Pi/Fst axes share this frame even with customized xy/width.
The existing gene y/height and each track's y/height settings remain independent.
Haplotype matrix columns remain equally spaced; connectors link their centers
to the true physical positions and the lower edge of the gene features.
No new parameter or output column is introduced. This replaces the former
inconsistent scaling of the gene body and the extended window region.

## Existing window tracks

The pi and Fst tracks display existing window CSV values; they do not estimate
statistics. Inputs from `hastat view -t pi` / `-t fst` can be used directly.
Required columns are `start,end,pi` or `start,end,fst`; `chrom` and `groups` are
preserved as literal strings (including numeric chromosome names and `NA`).
Use `chrom` for files containing multiple chromosomes. If it is absent, the file
is interpreted on the current gene's chromosome, retaining the earlier interface.

Coordinates are positive, 1-based inclusive integers. Windows must have
end >= start and be unique within chromosome/group. Overlapping sliding windows
are allowed. Rows fully contained in the selected gene region (including the
requested flanks) are displayed; partially overlapping windows are not cropped
or recomputed. Within each group, windows are sorted by midpoint.

The x-axis spans the selected region's outer base edges as defined above.
Each value has a point marker, so a single window or finite value isolated
by NaNs remains visible. Lines/fills retain NaN gaps. Missing values are never
converted to zero. Negative Fst values are retained. Pi is multiplied by 1000
only for the existing pi (10^-3) axis; the input CSV is unchanged. Infinite values,
negative pi, duplicate windows and invalid coordinates are rejected.

Input paths in a gene TOML are relative to that TOML's directory, or absolute.
An absent optional input setting omits that track; a supplied unreadable file is
an error. The separate `hap_color_file` remains a headerless haplotype-to-group
mapping, not the sample-to-haplotype file emitted by `view -t group`.

Scope: these repairs cover CSV ingress, display clustering, coordinate geometry
and window tracks. They do not change the existing unweighted
haplotype composition tracks. The latter are not sample-weighted population allele
frequencies. Native tool colors and user-defined layouts are retained. QA figures
are synthetic tool examples, not manuscript figure panels.
