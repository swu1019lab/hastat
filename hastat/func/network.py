# -*- coding: utf-8 -*-
# @Time    : 2024/12/22 10:00
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : network.py

from hastat.dataio.network import save_network
import time
import pandas as pd
from hastat.log.logger import logger
from hastat.utils.network import HaplotypeNetwork


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
    
    import csv
    with open(args.hap_table, encoding='utf-8-sig', newline='') as handle:
        columns = next(csv.reader(handle), [])
    if columns[:4] != ['chrom', 'pos', 'ref', 'alt'] or len(columns) < 5:
        raise ValueError("Haplotype table must start with chrom,pos,ref,alt and contain haplotype columns")
    if len(set(columns)) != len(columns):
        raise ValueError("Haplotype table column names must be unique")
    hap_data = pd.read_csv(args.hap_table, header=0)
    hap_data = hap_data.iloc[:, 4:]
    logger.info(f"Loaded haplotype table with {hap_data.shape[0]} loci and {hap_data.shape[1]} haplotypes")

    # Initialize network analysis
    network = HaplotypeNetwork()
    
    # Add edge data (calculate Hamming distances)
    network.set_haplotypes(hap_data)
    logger.info(f"Calculated Hamming distances for {hap_data.shape[1]} haplotypes")

    # Create network using specified method
    network.create_network(method=args.method)
    logger.info(f"Created haplotype network using {args.method} method")
    
    # Save network data if output file specified
    if args.out:
        save_network(network.network, args.out)

    end_time = time.time()
    logger.info("network function runs {:.2f} seconds".format(end_time - start_time))
