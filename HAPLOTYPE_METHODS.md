# Single-VCF genotype profiles and counts

`view -t table`, `group` and `freq` now share the same profile builder for
single regions and multiple loci. The builder is also used by `compare` after
its existing cross-VCF alignment and complete-locus selection. This patch is
deterministic grouping/counting; no statistical test, null distribution, P
value or new significance criterion is introduced.

## Profile definition

Apply the existing VCF reader and its requested `--het`, `--maf` and
`--max_missing` filters, evaluated over all VCF samples. Then exclude every
remaining locus with a missing diploid call in any input sample. Partial GTs
are whole missing calls under the existing reader policy. All samples use the
same complete loci; missing values do not create extra profiles. `geno` output
and its missing values remain unchanged.

The complete-locus requirement is stricter than a permissive `--max_missing`
setting. Allowing missing records through the reader does not make them usable
for grouping. Extensive missingness may leave no usable loci; then table and
group are empty, and freq has no assigned-sample counts. No universal Hap1 is
invented from an empty sequence. Retained identical profiles may differ at
discarded loci; full-sequence equivalence is not asserted.

These Hap labels describe unphased diploid 0/1/2 genotype vectors. They do not
resolve chromosome haplotypes or phase. Loci sort by chromosome string, position,
REF and ALT. Merge overlapping/adjacent requested intervals before extraction;
reject duplicate variant records instead of weighting them twice. Samples
retain their original VCF names. Rank profiles by abundance descending, then
genotype tuple descending, matching the previously repaired compare rule.
`table` columns and `group` labels always derive from that same mapping.

## Frequency means sample counts

The existing `freq` command produces **integer counts**, not proportions. This
meaning is retained. Ungrouped output has `haplotypes,samples` columns; the
sum of samples is the number of assigned VCF individuals. Grouped output has
`haplotypes,<group names>` with one flat header row, without the old redundant
outer `samples` header. Values are counts of individuals in each group/profile.
No percentage, denominator, status or diagnostic column is added.

The old multi-locus method grouped the two text columns of the assignment table
again instead of counting individuals, and ignored the provided group mapping.
It is replaced with the same counting routine as single-region freq.

Group metadata has two columns, sample and population. Preserve names such as
`001` and `NA` as strings. Every listed sample must exist in the VCF and appear
once, with a non-empty population name. Reject duplicate/unknown/empty entries
instead of multiplying or silently dropping samples. `haplotypes` is reserved
as a population label because it names the output row identifier.

A metadata file may intentionally list a subset of VCF samples. It restricts
counting only: profiles, complete-locus selection and Hap numbering still use
the entire VCF, keeping table/group/freq consistent. Unobserved haplotypes in
that subset are absent from freq; population columns are sorted and counts
retain the original global Hap IDs. With at least one usable locus, each column
sum equals the number of listed samples in that population. The subset scope
is recorded in the log, not additional result fields.

## Empty outputs and validation

Normal file types remain `.table.csv`, `.group.csv` and `.freq.csv`. Empty
results retain headers and replace stale results at an existing output prefix.
No extra result file or JSON sidecar is created.

`tests/test_haplotype.py` uses an independent Python dictionary/Counter oracle
and actual indexed VCF files. It verifies exact sample partitions, table/label
mapping and frequency totals, single/multi-region equivalence, grouping subsets,
filters, missingness, duplicates, literal IDs, ties, CLI serialization, repeat
runs and network input. Existing compare/stat/pi/D/Fst/network tests run in the
same PBS allocation. `validate_haplotype.py` records fixed synthetic before/after
cases under `qa/20261031_haplotype`; these are developer audit artifacts.

For genotype/phase encoding, see the [official VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf).
Complete-locus deletion and deterministic ranking are hastat software policies,
not requirements of the VCF specification. Numerical agreement validates the
implementation, not real-population phasing accuracy or large-cohort performance.
