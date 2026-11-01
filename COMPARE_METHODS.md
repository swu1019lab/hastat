# Cross-VCF genotype-profile comparison

This module performs deterministic alignment and grouping, not a statistical
test. It estimates no P values, ancestry, phasing, or population-specific allele
frequencies from missing records. The patch introduces no new statistical test.

## Alignment and available information

Input VCFs must already refer to the same reference assembly and chromosome
naming system, with consistently represented alleles. Align records by the
entire `(chrom,pos,ref,alt)` key, not POS alone. VCF sequence alleles are
case-insensitive, so nucleotide strings are uppercased for matching; symbolic
allele IDs remain case-sensitive. No strand conversion, REF/ALT swap, indel
normalization or liftover is inferred. Inconsistent SNP REF bases at a shared
chromosome/position cause an explicit error requiring input reconciliation.
Different ALT alleles with the same REF remain different records.

The shared reader retains biallelic diploid records. Missing or partial GT calls
are encoded as -1. A record absent from an input is not filled as 0/0. Its
absence may reflect filtering, coverage or VCF representation rather than a
biological absence. The existing `variant_class` fields `present_in_*`,
`present_in_all` and `population_specific` describe **record availability in
the supplied VCFs**, not proof that an allele is shared/private biologically.
An ALT allele need not be observed among the genotypes merely because the
VCF record exists.

Overlapping/adjacent requested intervals on each chromosome are merged before
fetching. Duplicate canonical variant records in one input cause an error;
overlapping requested regions do not duplicate a record. Genomic output order
is chromosome string, position, REF, ALT, independent of Python hash order.

## Common complete loci and grouping

Start with records present in every supplied VCF. Remove any of those loci
with a missing genotype in any included sample. All populations and the merged
analysis use exactly the same remaining loci. There is no imputation, no
missing-state haplotype, and no per-population selection of different sites.
The excluded count is logged. Extensive missingness may remove most/all loci;
the result describes only retained loci and can merge profiles distinguished
solely at excluded positions. This is not evidence that their full sequences
are identical.

Group individuals by their complete 0/1/2 dosage vectors over those loci.
These are **unphased multilocus genotype profiles**: existing Hap/MergedHap
names are retained as the tool's labels, not a claim of chromosome phase.
No complete shared loci means no assignments, not one universal haplotype.
Every input population remains represented in summaries, including its true
input sample count when it has zero variant records.

Rank profiles by sample count descending. Break equal-count ties by their
genotype tuple descending over the sorted loci. Local labels use `<Pop>HapN`;
merged labels use `MergedHapN`. Compare populations via the merged profiles:
`PopAHap1` and `PopBHap1` need not denote the same sequence. IDs are dataset-local
ranks and can change with sample composition or retained loci.

The established `_1`, `_2`, ... suffix on each sample is retained, identifying
the VCF's input order. Reordering VCF arguments changes those suffixes even
though the associated biological profile grouping is unchanged. Downstream
phenotype/group files must use the corresponding suffixed sample names.

## Validation and output behavior

Require at least two inputs, one unique non-empty population name per input,
and filename-safe labels. `merged` is reserved to prevent internal key collisions.
Non-default `--het`, `--maf` and `--max_missing` are rejected for `compare`;
they previously had no effect. Preprocess inputs explicitly if such filtering
is needed, remembering that filtering changes record availability.

Keep existing output types: population/merged `.table.csv` and `.group.csv`,
`.comparison_summary.csv`, and `.variant_class.csv`. Group files use the normal
two headers `samples,haplotypes` (fixing the accidental `index` name).
Empty haplotype results have headers and no assignments; they replace stale
results at the same prefix. No status columns, additional result types or JSON
sidecars are added. Existing summary rows and presence columns remain; `Core
variants` counts shared records **before** complete-locus removal, while the
haplotype tables show the actual retained loci. Internal validation records are
kept separately in `qa/20261031_compare/`.

An independent dictionary/Counter oracle verifies locus intersection, missing
call exclusion, individual partitions and count/rank mapping. Real indexed VCF
tests cover empty inputs, duplicates, ref/alt mismatches, sample/locus order,
hash-seed independence, CLI output and network loading. All earlier statistical
regressions run as well. Synthetic numerical/logic checks do not establish
biological equivalence or performance on large cohorts.

Format reference: [official VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf),
REF/ALT and genotype encoding/phase sections. The complete-locus rule and rank
tie breaker above are explicit hastat software policies, not requirements of
the VCF standard or newly proposed statistical estimators.
