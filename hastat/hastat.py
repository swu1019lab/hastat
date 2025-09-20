# -*- coding: utf-8 -*-
# @Time    : 2024/3/30 10:27
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : hastat.py


import argparse
import logging
from hastat.log.logger import logger
from .func import view, stat, plot, gwas, network


def main():
    parser = argparse.ArgumentParser(
        description='A package for gene haplotype analysis in natural populations',
        epilog="Designed on 02/22/2023 by Xiaodong Li"
    )

    parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')
    parser.add_argument("--log", type=str, help="The log file name (default: stdout)")
    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', help='additional help',
                                       dest='command')
    # Add subparsers for each command
    # view command
    parser_view = subparsers.add_parser('view', help='View and analyze the haplotypes data of genes or target regions')
    exclusive_group = parser_view.add_mutually_exclusive_group(required=True)
    parser_view.add_argument('-v', '--vcf', type=str, nargs='+', required=True, 
                            help='The VCF file(s) containing the genotype data. For compare type, multiple VCF files can be provided')
    parser_view.add_argument('-a', '--gff', type=str, help='The GFF file containing the gene annotation')
    exclusive_group.add_argument('-r', '--region', type=str,
                                 help='The region of the gene to be analyzed, format: chr:start-end')
    exclusive_group.add_argument('-i', '--gene_id', type=str,
                                 help='The gene ID for the target region, which should be provided with the GFF file')
    exclusive_group.add_argument('--homo', type=str,
                                 help='Multiple homologous gene IDs separated by comma, e.g., gene1,gene2,gene3')
    parser_view.add_argument('-u', '--upstream', type=int, default=0,
                             help='The upstream distance of the gene (default: 0)')
    parser_view.add_argument('-d', '--downstream', type=int, default=0,
                             help='The downstream distance of the gene (default: 0)')
    parser_view.add_argument('-t', '--type', type=str, choices=['geno', 'table', 'group', 'freq', 'pi', 'fst', 'compare'],
                             default='group', help='The data type to be analyzed (default: %(default)s)')
    parser_view.add_argument('-g', '--group', type=str,
                             help='A csv file containing the custom population groups of samples with two columns: sample name and group information')
    parser_view.add_argument('-w', '--size', type=int, default=1,
                             help='The size of the window for pi and fst analysis (default: %(default)s) (unit: bp)')
    parser_view.add_argument('-s', '--step', type=int, default=1,
                             help='The step size of the window for pi and fst analysis (default: %(default)s) (unit: bp)')
    parser_view.add_argument('-o', '--out', type=str, help='The output prefix for result files (default: stdout)')
    parser_view.add_argument('-n', '--pop_names', type=str, nargs='+',
                            help='Population names for compare analysis (must match the number of VCF files)')
    parser_view.add_argument('--het', type=float, default=None,
                            help='Heterozygosity rate threshold for filtering variants. Variants with heterozygosity rate greater than this value will be filtered out. 0=only homozygous sites, 1=only heterozygous sites (default: no filtering)')

    # stat command
    parser_stat = subparsers.add_parser('stat', help='Perform gene haplotype statistical test for related-traits')
    parser_stat.add_argument('-g', '--group', type=str, required=True,
                             help='A csv file containing the haplotype groups of samples with two columns: samples and haplotypes')
    parser_stat.add_argument('-p', '--pheno', type=str, required=True,
                             help='A csv file containing the phenotype data of samples with two columns at least: samples and trait name')
    parser_stat.add_argument('-s', '--min_hap_size', type=int, default=10,
                             help='The minimum sample size for each haplotype (default: %(default)s)')
    parser_stat.add_argument('-a', '--annotate', type=str, default='gene',
                             help='The annotation of haplotypes (default: %(default)s)')
    parser_stat.add_argument('-m', '--method', type=str, choices=['TukeyHSD', 'AllPairTest'], default='TukeyHSD',
                             help='The method for multiple comparisons (default: %(default)s)')
    parser_stat.add_argument('-o', '--out', type=str, help='The output csv file name (default: stdout)')

    # plot command with subcommands
    parser_plot = subparsers.add_parser('plot', help='Visualize gene haplotypes analysis results')
    plot_subparsers = parser_plot.add_subparsers(title='plot types', description='valid plot types', 
                                                help='different plot types', dest='plot_type')

    # bar plot subcommand
    parser_bar = plot_subparsers.add_parser('bar', help='Create bar plot for haplotype distribution')
    parser_bar.add_argument('--sample_hap', type=str, required=True,
                           help='CSV file with sample haplotype data (samples, haplotypes)')
    parser_bar.add_argument('--sample_group', type=str, required=True,
                           help='CSV file with sample group data (samples, group1, group2, ...)')
    parser_bar.add_argument('--group_index', type=int, default=1,
                           help='Index of group column to use (default: %(default)s)')
    parser_bar.add_argument('--calc_percentage', action='store_true',
                           help='Calculate percentage instead of count')
    parser_bar.add_argument('--haplotypes', type=str, nargs='+',
                           help='Haplotypes to be analyzed, if not provided, all haplotypes will be analyzed')
    parser_bar.add_argument('--x_label', type=str, default='Groups',
                           help='X-axis label (default: %(default)s)')
    parser_bar.add_argument('--y_label', type=str, default='Count',
                           help='Y-axis label (default: %(default)s)')
    parser_bar.add_argument('-o', '--output', type=str, required=True,
                           help='Output figure file path')
    parser_bar.add_argument('--width', type=float, default=5,
                           help='Figure width in inches (default: %(default)s)')
    parser_bar.add_argument('--height', type=float, default=4,
                           help='Figure height in inches (default: %(default)s)')

    # pie plot subcommand  
    parser_pie = plot_subparsers.add_parser('pie', help='Create pie plot for haplotype distribution')
    parser_pie.add_argument('--sample_hap', type=str, required=True,
                           help='CSV file with sample haplotype data (samples, haplotypes)')
    parser_pie.add_argument('--sample_group', type=str, required=True,
                           help='CSV file with sample group data (samples, group1, group2, ...)')
    parser_pie.add_argument('--group_index', type=int, default=1,
                           help='Index of group column to use (default: %(default)s)')
    parser_pie.add_argument('--calc_percentage', action='store_true',
                           help='Calculate percentage instead of count')
    parser_pie.add_argument('--haplotypes', type=str, nargs='+',
                           help='Haplotypes to be analyzed, if not provided, all haplotypes will be analyzed')
    parser_pie.add_argument('-o', '--output', type=str, required=True,
                           help='Output figure file path')
    parser_pie.add_argument('--width', type=float, default=5,
                           help='Figure width in inches (default: %(default)s)')
    parser_pie.add_argument('--height', type=float, default=4,
                           help='Figure height in inches (default: %(default)s)')

    # box plot subcommand
    parser_box = plot_subparsers.add_parser('box', help='Create box plot for haplotype-phenotype analysis')
    parser_box.add_argument('--sample_hap', type=str, required=True,
                           help='CSV file with sample haplotype data (samples, haplotypes)')
    parser_box.add_argument('--sample_phe', type=str, required=True,
                           help='CSV file with sample phenotype data (samples, trait1, trait2, ...)')
    parser_box.add_argument('--phe_index', type=int, default=1,
                           help='Index of phenotype column to use (default: %(default)s)')
    parser_box.add_argument('--comparisons', type=str, nargs='+', action='append',
                           help='Pairs of haplotypes for statistical comparison, e.g., --comparisons A B --comparisons A C')
    parser_box.add_argument('--method', type=str, choices=['t-test', 'mannwhitneyu', 'welch'], default='t-test',
                           help='Statistical test method for comparisons (default: %(default)s)')
    parser_box.add_argument('--step_size', type=float, nargs='+',
                           help='Step sizes for comparison brackets')
    parser_box.add_argument('--haplotypes', type=str, nargs='+',
                           help='Haplotypes to be analyzed, if not provided, all haplotypes will be analyzed')
    parser_box.add_argument('-o', '--output', type=str, required=True,
                           help='Output figure file path')
    parser_box.add_argument('--width', type=float, default=5,
                           help='Figure width in inches (default: %(default)s)')
    parser_box.add_argument('--height', type=float, default=4,
                           help='Figure height in inches (default: %(default)s)')
    
    # network plot subcommand
    parser_network = plot_subparsers.add_parser('network', help='Create network plot for haplotype network analysis')
    parser_network.add_argument('--file', type=str, required=True,
                               help='The network file (output from network command)')
    parser_network.add_argument('--sample_hap', type=str,
                               help='CSV file with sample haplotype data (samples, haplotypes)')
    parser_network.add_argument('--sample_group', type=str,
                               help='CSV file with sample group data (samples, group1, group2, ...)')
    parser_network.add_argument('--group_index', type=int, default=1,
                               help='Index of group column to use (default: %(default)s)')
    parser_network.add_argument('-o', '--output', type=str, required=True,
                           help='Output figure file path')
    parser_network.add_argument('--width', type=float, default=10,
                           help='Figure width in inches (default: %(default)s)')
    parser_network.add_argument('--height', type=float, default=8,
                           help='Figure height in inches (default: %(default)s)')
    parser_network.add_argument('--layout', type=str, default='spring',
                           choices=['spring', 'circular', 'kamada_kawai', 'fruchterman_reingold', 'shell', 'spectral', 'random'],
                           help='Network layout algorithm (default: %(default)s)')
    parser_network.add_argument('-k', '--k', type=float, default=None,
                           help='k: Optimal distance between nodes in spring layout, default is 1/sqrt(n) where n is the number of nodes. small k for more compact layout, large k for more spread out layout.')
    parser_network.add_argument('--show_node_label', action='store_true', default=False,
                           help='Show haplotype labels on nodes (default: %(default)s)')
    parser_network.add_argument('--node_label_size', type=float, default=10,
                           help='Node label font size (default: %(default)s)')
    parser_network.add_argument('--node_label_color', type=str, default='black',
                           help='Node label color (default: %(default)s)')
    parser_network.add_argument('--show_edge_label', action='store_true', default=False,
                           help='Show edge weight labels (default: %(default)s)')
    parser_network.add_argument('--edge_label_symbol', type=str, default=None,
                           help='Symbol to multiply with edge weights (e.g., "|" shows weight*symbol)')
    parser_network.add_argument('--edge_label_size', type=float, default=8,
                           help='Edge label font size (default: %(default)s)')
    parser_network.add_argument('--edge_label_color', type=str, default='black',
                           help='Edge label color (default: %(default)s)')
    parser_network.add_argument('--node_scale_factor', type=float, default=1.0,
                           help='Node size scaling factor (default: %(default)s)')
    parser_network.add_argument('--node_color', type=str, nargs='*',
                           help='Custom colors for population groups. Provide one or more colors in hex format (e.g., #FF0000 #00FF00 #0000FF)')
    parser_network.add_argument('--node_alpha', type=float, default=0.8,
                           help='Node color transparency (default: %(default)s)')
    parser_network.add_argument('--edge_color', type=str, default='gray',
                           help='Edge color (default: %(default)s)')
    parser_network.add_argument('--edge_width', type=float, default=1.0,
                           help='Edge width (default: %(default)s)')
    parser_network.add_argument('--edge_alpha', type=float, default=0.6,
                           help='Edge transparency (default: %(default)s)')
    parser_network.add_argument('--edge_style', type=str, default='line',
                           choices=['line', 'curve'],
                           help='Edge drawing style: line (straight) or curve (curved) (default: %(default)s)')
    parser_network.add_argument('--seed', type=int, default=None,
                           help='Random seed for reproducible layouts (default: %(default)s)')
    parser_network.add_argument('--title', type=str, default='Haplotype Network',
                           help='Plot title (default: %(default)s)')
    parser_network.add_argument('--legend_label_spacing', type=float, default=2,
                           help='The vertical space between the legend entries, in font-size units. (default: %(default)s)')
    parser_network.add_argument('--legend_handle_text_pad', type=float, default=1,
                           help='The pad between the legend handle and text, in font-size units. (default: %(default)s)')

    # gene plot subcommand
    parser_gene = plot_subparsers.add_parser('gene', help='Create gene structure plot with haplotype information')
    parser_gene.add_argument('--gff', type=str, required=True,
                            help='The GFF file containing the gene annotation')
    parser_gene.add_argument('--genes', type=str, nargs='+', required=True,
                            help='The gene IDs for the target region, which should be provided with the GFF file')
    parser_gene.add_argument('--toml', type=str, required=True,
                            help='The toml file containing the gene visualization configuration')
    parser_gene.add_argument('--upstream', type=int, default=0,
                            help='The upstream distance of the gene (default: %(default)s)')
    parser_gene.add_argument('--downstream', type=int, default=0,   
                            help='The downstream distance of the gene (default: %(default)s)')
    parser_gene.add_argument('-o', '--output', type=str, required=True,
                            help='Output figure file path')
    parser_gene.add_argument('--width', type=float, default=12.0,
                            help='Figure width in inches (default: %(default)s)')
    parser_gene.add_argument('--height', type=float, default=8.0,
                            help='Figure height in inches (default: %(default)s)')

    # gwas command
    parser_gwas = subparsers.add_parser('gwas', help='Perform GWAS analysis using GEMMA wrapper')
    parser_gwas.add_argument('-c', '--config', type=str, required=True, help='The configuration file for GWAS')

    # network command
    parser_network = subparsers.add_parser('network', help='Perform gene haplotype network analysis')
    parser_network.add_argument('-t', '--hap_table', type=str, required=True, 
                               help='A CSV file containing the haplotype table (output from view -t table)')
    parser_network.add_argument('-m', '--method', type=str, choices=['MST', 'MSN'], default='MST',
                               help='Network construction method: MST (Minimum Spanning Tree) or MSN (Minimum Spanning Network) (default: %(default)s)')
    parser_network.add_argument('-f', '--threshold_factor', type=float, default=1.2,
                               help='Threshold factor for MSN method, Edges with weight <= path_weight * threshold_factor will be added (default: %(default)s)')
    parser_network.add_argument('-o', '--out', type=str, default='network',
                               help='The output prefix for network analysis results (default: %(default)s)')

    # prepare the function specified by the subcommand
    parser_view.set_defaults(func=view.run)
    parser_stat.set_defaults(func=stat.run)
    parser_plot.set_defaults(func=plot.run)
    parser_bar.set_defaults(func=plot.run)
    parser_pie.set_defaults(func=plot.run)
    parser_box.set_defaults(func=plot.run)
    parser_gene.set_defaults(func=plot.run)
    parser_network.set_defaults(func=plot.run)
    parser_gwas.set_defaults(func=gwas.run)
    parser_network.set_defaults(func=network.run)

    args = parser.parse_args()

    # Set up logging
    if args.log:
        handler = logging.FileHandler(args.log, mode='w')
    else:
        handler = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # Call the function specified by the subcommand
    if args.func:
        args.func(args)


if __name__ == '__main__':
    main()
