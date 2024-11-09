# -*- coding: utf-8 -*-
# @Time    : 2024/9/30 18:07
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : plot.py

import time
import tomli
import matplotlib as mpl
mpl.use('Agg')

from hastat.log.logger import logger
from hastat.viz import bar, pie, box, network


viz_type = {
    'bar': bar.HapBar,
    'pie': pie.HapPie,
    'box': box.HapBox,
    'network': network.HapNetwork
}


def run(args):
    logger.info('Start to run plot function')
    start = time.time()

    with open(args.config, 'rb') as f:
        config = tomli.load(f)

    if viz_type.get(args.type) is None:
        logger.error("The plot type {} is not supported!!!".format(args.type))
    else:
        viz_type[args.type](config).plot()
    end = time.time()
    logger.info("plot function runs {:.2f} seconds\n".format(end - start))
