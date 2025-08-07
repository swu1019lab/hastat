# -*- coding: utf-8 -*-
# @Time    : 2024/9/30 18:09
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : gwas.py

import time
import tomli
from hastat.log.logger import logger
from hastat.utils.biotool import GemmaPipeline


def run(args):
    logger.info('Start to run gwas function')
    start = time.time()

    with open(args.config, 'rb') as f:
        config = tomli.load(f)

    gp = GemmaPipeline(gemma_path=config['tool']['gemma'], plink_path=config['tool']['plink'])
    logger.info(f"Gemma version: {gp.gemma.get_version()}")
    logger.info(f"Plink version: {gp.plink.get_version()}")
    gp.run(vcf_file=config['data']['vcf_file'],
           phe_file=config['data']['phe_file'],
           out_name=config['data']['out_name'],
           out_dir=config['data']['out_dir'],
           phe_num=config['data']['phe_num'])

    end = time.time()
    logger.info("gwas function runs {:.2f} seconds".format(end - start))
