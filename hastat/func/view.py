# -*- coding: utf-8 -*-
# @Time    : 2024/9/30 18:07
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : view.py

import time
import pandas as pd
from itertools import combinations
from hastat.log.logger import logger

# Set the display option to show all rows
pd.set_option('display.max_rows', None)

# Set the display option to show all columns
pd.set_option('display.max_columns', None)


def run(args):
    logger.info('Start to run view function')

    from hastat.utils import gene
    from hastat.dataio import gff, vcf

    start_time = time.time()
    gv = gene.GeneVariant(vcf.read(args.vcf))

    if args.region:
        chrom, start, end = args.region.replace(' ', '').replace('-', ':').split(':')
    elif args.gene_id:
        gf = gene.GeneFeature(gff.read(args.gff))
        chrom, start, end = gf.get_locus(args.gene_id)
    else:
        logger.error("The region or gene ID should be provided!!!")
        raise ValueError("The region or gene ID should be provided!!!")

    data = pd.DataFrame()
    if args.type == 'genotype':
        data = gv.get_geno_data(chrom, int(start), int(end))
    elif args.type == 'hap_table':
        data = gv.hap_table(chrom, int(start), int(end))
        data.reset_index(inplace=True)
    elif args.type == 'hap_group':
        data = gv.hap_groups(chrom, int(start), int(end))
    elif args.type == 'hap_freq':
        groups = pd.read_csv(args.group, header=0, names=['samples', 'groups']) if args.group else None
        data = gv.calc_hap_freq(chrom, int(start), int(end), groups)
        data.reset_index(inplace=True)
    elif args.type == 'pi':
        if not args.group:
            logger.error("The group file should be provided!!!")
            raise ValueError("The group file should be provided!!!")
        groups = pd.read_csv(args.group, header=0, names=['samples', 'groups'])
        data_list = []
        # list all groups to calculate Pi
        for group in groups.groups.unique():
            data = gv.get_pi_data(chrom, int(start), int(end),
                                  sample_list=groups.query('groups == @group').samples.tolist())
            data_list.append(pd.DataFrame(data, columns=['windows', 'n_bases', 'counts', 'pi']).assign(groups=group))
        data = pd.concat(data_list)
    elif args.type == 'fst':
        if not args.group:
            logger.error("The group file should be provided!!!")
            raise ValueError("The group file should be provided!!!")
        groups = pd.read_csv(args.group, header=0, names=['samples', 'groups'])
        data_list = []
        # list all pairwise comparisons to calculate Fst
        for pair in combinations(groups.groups.unique(), 2):
            data = gv.get_fst_data(chrom, int(start), int(end),
                                   groups.query('groups == @pair[0]').samples.tolist(),
                                   groups.query('groups == @pair[1]').samples.tolist())
            data_list.append(pd.DataFrame(data, columns=['windows', 'counts', 'fst']).assign(groups='_'.join(pair)))
        data = pd.concat(data_list)
    elif args.type == 'ld':
        data = gv.get_ld_data(chrom, int(start), int(end))
    else:
        logger.error("The data type {} is not supported!!!".format(args.type))

    if data.empty:
        logger.warning(f"The {args.type} data of target gene or region is empty!!!")

    if args.out:
        data.to_csv(args.out, index=False)
    else:
        print(data)

    end_time = time.time()
    logger.info("view function runs {:.2f} seconds\n".format(end_time - start_time))
