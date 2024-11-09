# -*- coding: utf-8 -*-
# @Time    : 2024/9/30 18:15
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : vcf.py

import os
from pysam import VariantFile
from hastat.hastat import logger


def read(file_path: str) -> VariantFile:
    """
    Load a vcf file contained genotype data of samples

    :param file_path: the path of the vcf file
    :return: a VariantFile object
    """
    logger.info("Reading the vcf file: %s", file_path)
    if os.path.exists(file_path):
        vcf_in = VariantFile(file_path)
        logger.info("The number of samples in the vcf file: %d", len(vcf_in.header.samples))
        return vcf_in
    else:
        logger.error("The file does not exist: %s", file_path)
        raise FileNotFoundError("The file does not exist: %s" % file_path)
