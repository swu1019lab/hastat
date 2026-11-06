#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"
mkdir -p ../images
hastat view -v gene.vcf.gz -r 1:1000001-3004000 -g populations.csv -t pi --size 50000 --step 50000 -o gene
hastat view -v gene.vcf.gz -r 1:1000001-3004000 -g populations.csv -t fst --size 50000 --step 50000 -o gene
hastat view -v gene.vcf.gz -r 1:1000001-3004000 -g populations.csv -t tajima_d --size 50000 --step 50000 -o gene
hastat plot gene --gff genes.gff3 --gene GeneA --up 2000 --down 2000 --width 15 --height 17.17300419287212 --toml gene_tracks.toml -o ../images/gene_tracks.pdf
hastat network -t gene.table.csv -m MST -o gene.mst.csv
hastat network -t gene.table.csv -m MSN -o gene.msn.csv
hastat plot network --file gene.mst.csv --sample_hap gene.group.csv --sample_group populations.csv --show_node_label --node_label_size 8.3 --show_edge_label --edge_label_size 8.4 --node_alpha 1 --width 6 --height 6 -o ../images/mst.png
hastat plot network --file gene.msn.csv --sample_hap gene.group.csv --sample_group populations.csv --show_node_label --node_label_size 8.3 --show_edge_label --edge_label_size 8.4 --node_alpha 1 --width 6 --height 6 -o ../images/msn.png
