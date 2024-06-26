# -*- coding: utf-8 -*-
# @Time    : 2024/3/30 10:27
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : hastat.py


import os
import argparse
import time
import logging
from multiprocessing import Pool, freeze_support, current_process

# Set up logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def worker(args, gid):
    """
    This function performs the main task of the script. It loads the data, retrieves the genotype of the gene,
    gets the haplotypes data of the gene, and performs statistical analysis on the haplotypes.

    :param args: An argparse.Namespace object containing the command-line arguments.
    :param gid: The id of the gene to be analyzed.
    :return: A string indicating the process name, the time taken, and the gene id.
    """
    try:
        from .dataset import DataSet
        from .stat import HapAnovaTest

        start = time.time()
        # Load data firstly
        db = DataSet(args.gen, args.ann, args.phe)
        # Get genotype object of gene
        geno = db.get_gene_geno(gid, upstream=args.up, downstream=args.down)
        if geno is None:
            logger.warning("\tFailed to retrieve results of {}!!!\n".format(gid))
        else:
            # Get haplotypes data of gene
            hap = geno.get_hap_data()
            hap.to_excel(file_path=os.path.join(args.dir, '{}.haplotypes.data.xlsx'.format(gid)))
            stat = HapAnovaTest(db.phe_in, hap.get_hap_groups(), min_hap_size=args.min_hap_size, annotate=gid)
            stat.export_data(csv_file=os.path.join(args.dir, '{}.haplotypes.stats.csv'.format(gid)))
            logger.info("\t{} is done!!!\n".format(gid))
        end = time.time()
        logger.info("{} runs {:.2f} seconds for {}\n".format(current_process().name, (end - start), gid))
    except Exception as e:
        logger.error(f"Error occurred: {e}")


def create_tasks(args):
    """
    This function creates a list of tasks to be processed.

    :param args: An argparse.Namespace object containing the command-line arguments.
    :return: A list of tasks.
    """
    TASKS = []
    with open(args.infile, 'rt') as fn:
        for line in fn:
            if line.strip("\n"):
                TASKS.append((args, line.strip("\n").split("\t")[0]))
    return TASKS


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
    parser.add_argument("--log", type=str, default='hastat.log', help="The log file (default: %(default)s)")
    parser.add_argument("--gen", required=True, type=str, help="A genotype file with VCF format")
    parser.add_argument("--ann", required=True, type=str, help="A gene annotation file with GFF format")
    parser.add_argument("--phe", required=True, type=str,
                        help="A phenotype file with CSV format (the first column is samples and others are phenotypes)")
    parser.add_argument("--process", type=int, default=1,
                        help="The number of processes to run (default: %(default)s)")
    parser.add_argument("--min-hap-size", dest="min_hap_size", type=int, default=10,
                        help="the minimum sample size for each haplotype, default is 10 (default: %(default)s)")
    parser.add_argument("--gene-up-stream", dest="up", type=int, default=1500,
                        help="up stream length of gene (default: %(default)s)")
    parser.add_argument("--gene-down-stream", dest="down", type=int, default=0,
                        help="down stream length of gene (default: %(default)s)")
    parser.add_argument('--out-dir', dest="dir", default='hap_results',
                        help="The directory of output file. (default: %(default)s)")
    parser.add_argument("infile", help="A gene list file with one gene id per line")

    args = parser.parse_args()

    # Set up logging
    handler = logging.FileHandler(args.log, mode='w')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # Create output directory
    if not os.path.exists(args.dir):
        os.makedirs(args.dir)

    # Create tasks
    TASKS = create_tasks(args)

    # Start the timer and running the tasks
    NUMBER_OF_PROCESSES = args.process
    logger.info('Creating pool with {} processes\n'.format(NUMBER_OF_PROCESSES))

    start = time.time()
    # Create Pool
    with Pool(NUMBER_OF_PROCESSES) as pool:
        try:
            pool.starmap(worker, TASKS)
        except Exception as e:
            logger.error(f"Error occurred during multiprocessing: {e}")
    end = time.time()

    logger.info("All subprocesses done within {:.2f} seconds.".format(end - start))


if __name__ == '__main__':
    freeze_support()
    main()
