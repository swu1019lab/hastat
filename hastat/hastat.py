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
from .dataset import DataSet
from .stat import HapAnovaTest

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def worker(args, out_dir, gid):
    """
    This function performs the main task of the script. It loads the data, retrieves the genotype of the gene,
    gets the haplotypes data of the gene, and performs statistical analysis on the haplotypes.

    :param args: An argparse.Namespace object containing the command-line arguments.
    :param out_dir: The directory where the output files will be saved.
    :param gid: The id of the gene to be analyzed.
    :return: A string indicating the process name, the time taken, and the gene id.
    """
    try:
        start = time.time()
        # Load data firstly
        db = DataSet(args.gen, args.ann, args.phe)
        # Get genotype object of gene
        geno = db.get_gene_geno(gid, upstream=args.up, downstream=args.down)
        if geno is None:
            logging.warning("\tFailed to retrieve results of {}!!!\n".format(gid))
        else:
            # Get haplotypes data of gene
            hap = geno.get_hap_data()
            hap.to_excel(file_path=os.path.join(out_dir, '{}.haplotypes.data.xlsx'.format(gid)))
            stat = HapAnovaTest(db.phe_in, hap.get_hap_ngroup())
            stat.to_excel(file_path=os.path.join(out_dir, '{}.haplotypes.stats.xlsx'.format(gid)))

        end = time.time()
        return "{} runs {:.2f} seconds for {}\n".format(current_process().name, (end - start), gid)
    except Exception as e:
        logging.error(f"Error occurred: {e}")


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
                TASKS.append((args, args.dir, line.strip("\n").split("\t")[0]))
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
    parser.add_argument('--version', action='version', version='%(prog)s 0.0.7')
    parser.add_argument("--gen", required=True, type=str, help="A genotype file with VCF format")
    parser.add_argument("--ann", required=True, type=str, help="A gene annotation file with GFF format")
    parser.add_argument("--phe", required=True, type=str,
                        help="A phenotype file with CSV format (the first column is samples and others are phenotypes)")
    parser.add_argument("--process", type=int, default=1, help="The number of processes to run (default: %(default)s)")
    parser.add_argument("--hap-min-num", type=int, default=0, help="The minimum number of haplotypes for analysis")
    parser.add_argument("--gene-up-stream", dest="up", type=int, default=1500, help="up stream length of gene")
    parser.add_argument("--gene-down-stream", dest="down", type=int, default=0, help="down stream length of gene")
    parser.add_argument('--out-dir', dest="dir", default='data',
                        help="The directory of output file. (default: %(default)s)")

    parser.add_argument("infile", help="A text file containing gene id per line only")

    args = parser.parse_args()

    out_dir = os.path.join(args.dir, time.strftime("%Y%m%d", time.localtime()))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Create tasks
    TASKS = create_tasks(args)

    # Start the timer and running the tasks
    start = time.time()
    NUMBER_OF_PROCESSES = args.process
    logging.info('Creating pool with {} processes\n'.format(NUMBER_OF_PROCESSES))
    # Create Pool
    with Pool(NUMBER_OF_PROCESSES) as pool:
        results = [pool.apply_async(worker, task) for task in TASKS]
        for r in results:
            logging.info(r.get())
    end = time.time()

    logging.info("All subprocesses done within {:.2f} seconds.".format(end - start))


if __name__ == '__main__':
    freeze_support()
    main()
