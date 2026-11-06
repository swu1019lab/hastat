# Complete graphical examples

The [gallery](gallery.html) shows the GeneA plot with 1 Mb flanks, haplotype
networks, and bar, pie and box plots.
The haplotype matrix contains 48 variants from **160 synthetic samples and eight
haplotypes**. The sample metadata identifies 49 Wild, 61 Landrace and 50 Cultivar
samples. Hap1–Hap8 contain 36/30/24/20/18/14/10/8 samples, respectively.
Each haplotype contains all three material types. The window estimates use
the example VCF and these same sample assignments.

After installing hastat, run from the repository root:

    bash docs/examples/commands.sh

The script creates `docs/images/gene_tracks.pdf` and separate `docs/images/mst.png`
and `docs/images/msn.png` network plots. Each network plotting call uses one file.
The gallery displays PNG previews and their plotting commands.
Run the tutorial commands below from `docs/examples`.

All ordinary `hastat plot` exports also embed TrueType fonts in PDF
(`pdf.fonttype=42`) and preserve text elements in SVG (`svg.fonttype=none`).
No additional option is required. For text editing in Adobe Illustrator, install
the fonts used by the plot or substitute available fonts. The supplied examples
use DejaVu Sans, included with Matplotlib. Equations and rotated labels may
import as multiple text objects.

## Inputs

| File | Contents |
|:---|:---|
| [gene.table.csv](examples/gene.table.csv) | chrom,pos,ref,alt followed by eight haplotype genotype-state columns; 48 rows. |
| [gene.group.csv](examples/gene.group.csv) | Unique samples,haplotypes assignments; supplies sample counts. |
| [populations.csv](examples/populations.csv) | Sample ID and material type (Wild, Landrace or Cultivar). Used for network pies, allele frequencies, window statistics, and bar and pie counts. |
| [phenotypes.csv](examples/phenotypes.csv) | samples,Trait numerical observations. |
| [genes.gff3](examples/genes.gff3) | GeneA, chromosome 1, bases 2,000,001–2,004,000, plus strand and five CDS blocks. |
| [hap_colors.csv](examples/hap_colors.csv) | Haplotype row colors: group1 contains Hap1–3, group2 contains Hap6–7, and group3 contains Hap4, Hap5 and Hap8. The CSV has no header. |
| [annotations.csv](examples/annotations.csv) | Positions and categories for the 48 plotted variants, with optional labels on two CDS sites. |
| [gene.pi.csv](examples/gene.pi.csv), [gene.fst.csv](examples/gene.fst.csv) | 41 windows from 1,000,001 to 3,004,000 per series: π for each material type and Fst for each pair of types. Window size and step are both 50 kb. |
| [gene.vcf](examples/gene.vcf) | Diploid genotypes for all 160 samples, including the 48 displayed variants and additional flanking sites used in the statistics tracks. |
| [gene.vcf.gz](examples/gene.vcf.gz) and [its tabix index](examples/gene.vcf.gz.tbi) | Indexed copy used directly by the gene figure. |
| [gene.tajima_d.csv](examples/gene.tajima_d.csv) | Tajima's D for the three material types, using the same 41 windows as π and Fst. The last window is clipped at 3,004,000. |
| [gene_tracks.toml](examples/gene_tracks.toml) | Complete gene-plot configuration, including frequency and window tracks. |
| [gene.mst.csv](examples/gene.mst.csv), [gene.msn.csv](examples/gene.msn.csv) | MST and MSN edges generated from the same 48-variant table. |

Coordinates are 1-based inclusive. The 48-site haplotype matrix uses homozygous states
0 (reference) and 2 (alternate); the tool also recognizes 1 (heterozygous) and -1
(missing). Synthetic sample identifiers are demonstration labels.

## 1. Population sample counts

![Population sample counts](images/bar.png)

```bash
hastat plot bar --sample_hap gene.group.csv --sample_group populations.csv --x_label Material --y_label Samples --width 10 --height 4 -o ../images/bar.pdf
```

Stack segments are counts of matched samples assigned to each haplotype.
For relative composition, add --calc_percentage and use --y_label Fraction:
that option returns values between 0 and 1.
[PDF](images/bar.pdf) · [SVG](images/bar.svg)

## 2. Population composition

![Population composition](images/pie.png)

```bash
hastat plot pie --sample_hap gene.group.csv --sample_group populations.csv --width 11 --height 4 -o ../images/pie.pdf
```

Each population has its own pie, with the matched sample total below it.
Haplotype colors are shared with the bar plot.
[PDF](images/pie.pdf) · [SVG](images/pie.svg)

## 3. Trait distributions and comparisons

![Trait distributions](images/box.png)

```bash
hastat plot box --sample_hap gene.group.csv --sample_phe phenotypes.csv --method TukeyHSD --comparisons Hap1 Hap2 --comparisons Hap1 Hap4 --comparisons Hap1 Hap8 --width 10 --height 5 -o ../images/box.pdf
```

Points show observations; boxes show medians and interquartile ranges.
Whiskers extend to the most extreme observations within 1.5 interquartile ranges.
The three displayed comparisons use Tukey-adjusted P values from the full
eight-haplotype family (28 pairwise comparisons), even though only three
brackets are requested. Symbols: *** < 0.001, ** < 0.01, * < 0.05, ns ≥ 0.05.
TukeyHSD assumes independent observations, approximately normal within-group
errors and common within-group variance.
[PDF](images/box.pdf) · [SVG](images/box.svg)

## 4. Minimum spanning tree and network

![MST and MSN at the same node positions](images/mst_msn.png)

The left panel shows one minimum spanning tree (seven edges). The right panel
adds three equal-cost alternatives: Hap1–Hap3 (16 differing sites), Hap1–Hap7
(16) and Hap6–Hap8 (18). Blue dashed lines distinguish those edges from the
displayed tree. The two panels share the same weighted Kamada–Kawai coordinates
and node size scale. Pie sectors show the fractions of Wild, Landrace and
Cultivar samples assigned to each haplotype. For example, Hap1 has 14 Wild,
14 Landrace and 8 Cultivar samples. Circle area is proportional to sample count,
using the nested 8-, 20- and 36-sample circles as a reference. Only samples
with a value in the selected metadata column contribute to the counts.
The network represents genotype differences, not ancestry or evolutionary time.
[PDF](images/mst_msn.pdf) · [SVG](images/mst_msn.svg)

## 5. Gene haplotypes with 1 Mb flanks

![Complete gene haplotype figure](images/gene_tracks.png)

Use the [recommended GeneA TOML](examples/gene_tracks.toml) with the supplied
[example inputs](examples/), including the indexed VCF. From that directory, run:

```bash
hastat plot gene --gff genes.gff3 --gene GeneA --up 2000 --down 2000 --width 15 --height 17.17300419287212 --toml gene_tracks.toml -o gene_tracks.pdf
```

The command selects GeneA, its GFF, both flanks and the canvas dimensions in inches.
`[default]` sets shared positions and heights within each gene's `bbox` (0 at the
bottom, 1 at the top). `[GeneA]` inherits those settings and supplies its own inputs.
Select a different gene from the same TOML by passing its ID to `--gene`.
The panel `bbox` sets the left and right edges of every track. `gene_y` sets the vertical center of the gene.
Each statistics track shows all windows for GeneA's chromosome in its own CSV.
The supplied windows extend 1 Mb on either side of GeneA. The gene structure
displays 2 kb of flanking sequence on each side.
With `layout = "fixed"`, increasing `pi_y` for GeneA moves only its Pi axis upwards.
`ld_y` anchors the top of the LD triangle. The VCF supplies data for AF and LD;
changing input paths does not select a layout. The panel order from top to bottom is:

1. Pi and Fst for the three populations and their pairwise comparisons, using
   the supplied 50 kb window estimates. Pi is displayed in units of 10⁻³.
   With `shade_gene = true`, a gray rectangle marks GeneA (2,000,001–2,004,000) on all
   three statistics tracks. Its flanking regions lie outside the rectangle.
2. Tajima's D for the same three material types and 50 kb windows. Points mark window
   midpoints and the dashed line marks zero.
3. The five dark-blue CDS blocks of a 4 kb gene with 2 kb flanks. A single
   `>` or `<` at the strand-direction end marks orientation, without a connecting
   tail outside the gene. Short labels identify selected variant positions.
4. An annotation strip above the 48 variant columns and eight clustered
   haplotype rows, with sample counts and haplotype category colors.
5. REF/ALT allele frequencies in Wild, Landrace and Cultivar. Pale blue denotes
   REF and dark blue denotes ALT. At 1:1,998,150, Hap2 and Hap3 carry ALT. The called
   ALT fractions are 40/98 (40.8%) in Wild, 46/122 (37.7%) in Landrace and
   22/100 (22.0%) in Cultivar.
6. Whole-population LD for the 48 displayed sites in the same column order.
   Each diamond shows pairwise dosage r² from the 160 VCF samples, with gray
   for undefined pairs. The triangle is automatically half as high as it is wide.

The file prefix requires gene.table.csv and discovers optional gene.group.csv.
The annotation strip distinguishes upstream, CDS, intron and downstream sites.
Its legend occupies one row directly above the gene and its callouts.
The optional label field highlights two positions. Short labels are placed in
up to two rows, and annotations at the same variant share a cell.
Set `ld_show = false` to hide LD while retaining VCF-derived AF. Use
`ld_cmap` and `ld_colorbar` to choose the palette and toggle the color scale.
When the group file is absent, labels omit sample counts. hap_color_file
controls the row colors. The legend below the tree uses the category names in
that file. All paths in the TOML are relative to
that TOML file.

Clustering uses average linkage on Hamming distances over common complete loci.
It describes profile similarity, not a phylogeny or a significance test.
The 4 kb gene occupies 50.0% of the 8 kb gene view and about 0.2% of each
2.004 Mb statistics track. The matrix columns stay equally spaced.

The example script reads the supplied sample-to-material table and runs `hastat view` for
Pi, Fst and Tajima's D on the same samples and genomic interval. Undefined
estimates appear as gaps. The gene plot then reads the three CSV files.
The Tajima's D CSV uses
`chrom,start,end,counts,tajima_d,groups`. In the linked TOML, `tajima_d_y`
and `tajima_d_height` position its axis beneath Pi and Fst.
[PDF](images/gene_tracks.pdf) · [SVG](images/gene_tracks.svg)


