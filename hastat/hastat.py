# -*- coding: utf-8 -*-
# @Time    : 2024/3/30 10:27
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : hastat.py


import argparse
import logging
from hastat.log.logger import logger
from .func import view, stat, plot, gwas


def main():
    """
    This function handles command-line arguments and calls the worker function. It creates a pool of processes
    and applies the worker function to each task in the TASKS list.

    :return: None
    """
    parser = argparse.ArgumentParser(
        description='A package for gene haplotype analysis',
        epilog="Designed on 02/22/2023 by Xiao dong Li"
    )

    parser.add_argument('--version', action='version', version='%(prog)s 0.0.4')
    parser.add_argument("--log", type=str, help="The log file name (default: stdout)")
    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', help='additional help',
                                       dest='command')
    # Add subparsers for each command
    # view command
    parser_view = subparsers.add_parser('view', help='View the haplotypes data of a gene or target region')
    exclusive_group = parser_view.add_mutually_exclusive_group(required=True)
    parser_view.add_argument('-v', '--vcf', type=str, required=True, help='The VCF file containing the genotype data')
    parser_view.add_argument('-a', '--gff', type=str, help='The GFF file containing the gene annotation')
    exclusive_group.add_argument('-r', '--region', type=str,
                                 help='The region of the gene to be analyzed, format: chr:start-end')
    exclusive_group.add_argument('-i', '--gene_id', type=str,
                                 help='The gene ID for the target region, which should be provided with the GFF file')
    parser_view.add_argument('-t', '--type', type=str, choices=['genotype', 'hap_table', 'hap_group', 'hap_freq'],
                             default='hap_group', help='The data type to be analyzed (default: %(default)s)')
    parser_view.add_argument('-g', '--group', type=str,
                             help='A csv file containing the custom groups of samples if the data type is hap_group')
    parser_view.add_argument('-o', '--out', type=str, help='The output csv file name (default: stdout)')

    # stat command
    parser_stat = subparsers.add_parser('stat', help='Perform statistical analysis on the haplotypes of a gene')
    parser_stat.add_argument('-g', '--group', type=str, required=True,
                             help='A csv file containing the haplotype groups of samples')
    parser_stat.add_argument('-p', '--pheno', type=str, required=True,
                             help='A csv file containing the phenotype data of samples')
    parser_stat.add_argument('-s', '--min_hap_size', type=int, default=10,
                             help='The minimum sample size for each haplotype (default: %(default)s)')
    parser_stat.add_argument('-a', '--annotate', type=str, default='gene',
                             help='The annotation of haplotypes (default: %(default)s)')
    parser_stat.add_argument('-m', '--method', type=str, choices=['TukeyHSD', 'AllPairTest'], default='TukeyHSD',
                             help='The method for multiple comparisons (default: %(default)s)')
    parser_stat.add_argument('-o', '--out', type=str, help='The output csv file name (default: stdout)')

    # plot command
    parser_plot = subparsers.add_parser('plot', help='Plot the haplotypes data of a gene')
    parser_plot.add_argument('-t', '--type', type=str, choices=['bar', 'pie', 'box', 'network'], required=True,
                             help='The plot type')
    parser_plot.add_argument('-c', '--config', type=str, required=True, help='The configuration file for plotting')

    # gwas command
    parser_gwas = subparsers.add_parser('gwas', help='Perform GWAS analysis using GEMMA/EMAX wrapper')
    parser_gwas.add_argument('-c', '--config', type=str, required=True, help='The configuration file for GWAS')

    # prepare the function specified by the subcommand
    parser_view.set_defaults(func=view.run)
    parser_stat.set_defaults(func=stat.run)
    parser_plot.set_defaults(func=plot.run)
    parser_gwas.set_defaults(func=gwas.run)

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
