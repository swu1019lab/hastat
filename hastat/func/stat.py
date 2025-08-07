# -*- coding: utf-8 -*-
# @Time    : 2024/9/30 18:08
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : stat.py

import time
import pandas as pd
from prettytable.colortable import ColorTable, Themes
from hastat.log.logger import logger


def run(args):
    logger.info('Start to run stat function')
    start = time.time()

    from hastat.utils import stat
    pheno = pd.read_csv(args.pheno)
    groups = pd.read_csv(args.group, header=0, names=['samples', 'haplotypes'])
    res = stat.HapAnovaTest(pheno, groups, args.min_hap_size, args.annotate, args.method)
    if args.out:
        pd.concat(res.anova_data['comparisons']).to_csv(args.out, index=False)
    else:
        df = pd.concat(res.anova_data['comparisons'])
        table = ColorTable(theme=Themes.OCEAN)
        table.field_names = df.columns.tolist()
        for row in df.itertuples(index=False):
            table.add_row(row)
        print(table)
    end = time.time()
    logger.info("stat function runs {:.2f} seconds".format(end - start))
