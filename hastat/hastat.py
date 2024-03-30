# -*- coding: utf-8 -*-
# @Time    : 2024/3/30 10:27
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : hastat.py


import os
import argparse
import time
from multiprocessing import Pool, freeze_support, current_process
from .dataset import DataSet
from .stat import HapAnovaTest


def worker(args, out_dir, gid):
    start = time.time()
    # Load data firstly
    db = DataSet(args.gen, args.ann, args.phe)
    # Get genotype object of gene
    geno = db.get_gene_geno(gid, upstream=args.up, downstream=args.down)
    if geno is None:
        print("\tFailed to retrieve results of {}!!!\n".format(gid))
    else:
        # Get haplotypes data of gene
        hap = geno.get_hap_data()
        hap.to_excel(file_path=os.path.join(out_dir, '{}.haplotypes.data.xlsx'.format(gid)))
        stat = HapAnovaTest(db.phe_in, hap.get_hap_ngroup())
        stat.to_excel(file_path=os.path.join(out_dir, '{}.haplotypes.stats.xlsx'.format(gid)))

    end = time.time()
    return "{} runs {:.2f} seconds for {}\n".format(current_process().name, (end - start), gid)


def main():
    parser = argparse.ArgumentParser(
        description='A package for gene haplotype analysis',
        epilog="Designed on 02/22/2023 by Xiao dong Li"
    )
    parser.add_argument('--version', action='version', version='%(prog)s 0.0.7')
    parser.add_argument("--log", type=str, default="hap_stat.log", help="Output log file (default: %(default)s)")
    parser.add_argument("--gen", required=True, type=str, help="A genotype file with VCF format")
    parser.add_argument("--ann", required=True, type=str, help="A gene annotation file with GFF format")
    parser.add_argument("--phe", required=True, type=str,
                        help="A phenotype file with CSV format (the first column is samples and others are phenotypes)")
    parser.add_argument("--process", type=int, default=4, help="")
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
    TASKS = []
    with open(args.infile, 'rt') as fn:
        for line in fn:
            if line.strip("\n"):
                TASKS.append((args, out_dir, line.strip("\n").split("\t")[0]))

    start = time.time()
    NUMBER_OF_PROCESSES = args.process
    print('Creating pool with {} processes\n'.format(NUMBER_OF_PROCESSES))
    # Create Pool
    with Pool(NUMBER_OF_PROCESSES) as pool:
        results = [pool.apply_async(worker, task) for task in TASKS]
        for r in results:
            print(r.get())
    end = time.time()

    print("All subprocesses done within {:.2f} seconds.".format(end - start))


if __name__ == '__main__':
    freeze_support()
    main()
