# -*- coding: utf-8 -*-
# @Time    : 2024/9/30 18:07
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : plot.py

import time
import tomli
from hastat.log.logger import logger
from hastat.viz import bar, pie, box, network, gene
from hastat.manage import Manager
import matplotlib as mpl
mpl.use('Agg')

# Initialize the manager and register plot types
manager = Manager()
manager.register('bar', bar.HapBar)
manager.register('pie', pie.HapPie)
manager.register('box', box.HapBox)
manager.register('network', network.HapNetwork)
manager.register('gene', gene.HapGene)


def run(args):
    logger.info('Start to run plot function')
    start = time.time()

    with open(args.config, 'rb') as f:
        config = tomli.load(f)

    try:
        plot_instance = manager.use(args.type, config)
        plot_instance.plot()
    except ValueError as e:
        logger.error(e)

    end = time.time()
    logger.info("plot function runs {:.2f} seconds\n".format(end - start))
