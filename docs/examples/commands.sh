#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"
mkdir -p ../images
python -c "import pysam; pysam.tabix_compress('gene.vcf', 'gene.vcf.gz', force=True); pysam.tabix_index('gene.vcf.gz', preset='vcf', force=True)"
hastat view -v gene.vcf.gz -r 1:501-8500 -g populations.csv -t tajima_d --size 1000 --step 500 -o gene
hastat network -t gene.table.csv -m MSN -o gene.network.csv
hastat plot bar --sample_hap gene.group.csv --sample_group populations.csv --x_label Population --y_label Samples --width 10 --height 4 -o ../images/bar.pdf
hastat plot pie --sample_hap gene.group.csv --sample_group populations.csv --width 11 --height 4 -o ../images/pie.pdf
hastat plot box --sample_hap gene.group.csv --sample_phe phenotypes.csv --method TukeyHSD --comparisons Hap1 Hap2 --comparisons Hap1 Hap4 --comparisons Hap1 Hap8 --width 10 --height 5 -o ../images/box.pdf
hastat plot network --file gene.network.csv --sample_hap gene.group.csv --sample_group populations.csv --show_node_label --show_edge_label --layout spring --seed 31 --node_scale_factor 550 --node_alpha 1 --width 14 --height 9 --title 'Haplotype network' -o ../images/network.pdf
hastat plot gene --gff genes.gff3 --genes GeneA --toml gene_tracks.toml --upstream 2000 --downstream 2000 --width 15 --height 13.2 -o ../images/gene_tracks.pdf
