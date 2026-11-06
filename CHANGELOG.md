# Changelog

## hastat 2.0.0

- Show effective sample counts and overall test P values on box plots; add --x_rotation for label angles.
- Add hap_min_size to filter displayed haplotypes by sample count while retaining the full population for AF and LD.
- Automatically draw CDS and UTR blocks, falling back to exons when CDS annotations are absent; remove the gene_feature setting.
- Place terminal strand markers closer to the gene boundary.
- Refresh the box and gene example figures and their documentation.

- Reorganize Python APIs for variant data, haplotypes, population statistics,
  and haplotype-trait association.
- Organize gene visualization into gene.py and tracks.py.
- Add population-aware REF/ALT frequency tracks based on called allele
  counts, and integrate whole-population LD heatmaps into gene plots.
- Add compact variant annotation bands with selected labels.
- Refine shared and per-gene TOML settings, align tracks to bbox, and
  support independent window-statistic ranges with optional gene shading.
- Use --gene for single-gene visualization and --up/--down for flanks.
- Improve network plots with weighted Kamada–Kawai layouts, population
  pie nodes, proportional node areas, and nested sample-count legends.
- Refine gene colors, single terminal strand markers, statistical labels,
  and pie-chart typography.
- Update complete examples, the recommended TOML configuration,
  and the PNG gallery.

- Add searchable web documentation with installation instructions, CLI and TOML references, detailed examples, and a PNG gallery.
- Simplify README examples and link to the web documentation.
- Render network genotype differences as perpendicular ticks on straight and curved edges.
- Use Arial as the default plotting font and refresh all example figures.
- Remove the GEMMA-based GWAS command and its bundled configuration.

## 1.1.0 (2026-11-01)

- Corrected ANOVA and multiple comparisons, including cases with small groups, zero variance or inestimable contrasts. Use `WelchBonferroni` in place of `AllPairTest` for Welch’s ANOVA followed by pairwise Welch t-tests with Bonferroni correction.
- Corrected calculations and window boundaries for nucleotide diversity (π), Tajima’s D and Fst. Undefined estimates are reported as `NaN`.
- Corrected haplotype grouping across regions, sample counts, variant matching across files and minimum spanning network construction.
- Gene haplotype plots accept a file prefix and show sample counts when the matching `.group.csv` exists. Clustering uses Hamming distance and average linkage by default. Gene plots show strand direction, flanking coordinates and window statistics. Use `hap_color_file` in place of `hap_group_file` for haplotype category colors.
- Network nodes cover their edges. Population colors match the other plots, and the size legend shows sample counts.
- PDF exports embed TrueType fonts and SVG exports preserve text. The gallery includes complete figures and the commands to reproduce them.
