# -*- coding: utf-8 -*-
# @Time    : 2024/12/22 10:00
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : network.py

import time
import pandas as pd
from hastat.log.logger import logger
from hastat.utils.gene import GeneHapNetwork


def run(args):
    """
    Run haplotype network analysis

    :param args: command line arguments
    :return: None
    """
    logger.info('Start to run network function')
    start_time = time.time()

    # Read haplotype table data
    if not args.hap_table:
        logger.error("Haplotype table file should be provided!")
        raise ValueError("Haplotype table file should be provided!")
    
    hap_data = pd.read_csv(args.hap_table, header=0)
    hap_data = hap_data.iloc[:, 4:] # exclude the first 4 columns: chrom, pos, ref, alt
    logger.info(f"Loaded haplotype table with {hap_data.shape[0]} loci and {hap_data.shape[1]} haplotypes")

    # Initialize network analysis
    network = GeneHapNetwork()
    
    # Add edge data (calculate Hamming distances)
    network.add_edge(hap_data)
    logger.info(f"Calculated Hamming distances for {hap_data.shape[1]} haplotypes")

    # Create network using specified method
    if args.method == 'MST':
        network.create_network(method=args.method)
    elif args.method == 'MSN':
        network.create_network(method=args.method, threshold_factor=args.threshold_factor)
    logger.info(f"Created haplotype network using {args.method} method")
    
    # Save network data if output file specified
    if args.out:
        network.save_network(args.out)

    end_time = time.time()
    logger.info("network function runs {:.2f} seconds".format(end_time - start_time))
