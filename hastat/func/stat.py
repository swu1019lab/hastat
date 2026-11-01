# -*- coding: utf-8 -*-
# @Time    : 2024/9/30 18:08
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : stat.py

import time
import json
from hastat.utils.sample_data import read_sample_csv, read_assignments
from prettytable.colortable import ColorTable, Themes
from hastat.log.logger import logger


def run(args):
    logger.info('Start to run stat function')
    start = time.time()

    from hastat.utils import stat
    pheno = read_sample_csv(args.pheno)
    groups = read_assignments(args.group, label='haplotypes')
    res = stat.HapAnovaTest(pheno, groups, args.min_hap_size, args.annotate, args.method,
                          alpha=getattr(args, 'alpha', 0.05))
    if args.out:
        res.comparisons.to_csv(args.out, index=False, float_format='%.17g')
    else:
        df = res.comparisons
        table = ColorTable(theme=Themes.OCEAN)
        table.field_names = df.columns.tolist()
        for row in df.itertuples(index=False):
            table.add_row(row)
        print(table)
        for result in res.anova_data['omnibus']:
            print(json.dumps(result, allow_nan=False))
    end = time.time()
    logger.info("stat function runs {:.2f} seconds".format(end - start))
