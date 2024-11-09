# -*- coding: utf-8 -*-
# @Time    : 2024/9/30 18:14
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : gff.py

import os
import gffutils
from hastat.hastat import logger


def read(file_path: str) -> gffutils.FeatureDB:
    """
    Load a gff file contained gene annotation data of samples

    :param file_path: the path of the gff file
    :return: a FeatureDB object
    """
    logger.info("Reading the gff file: %s", file_path)
    try:
        if os.path.exists(file_path + '.sqlite3'):
            gff_in = gffutils.FeatureDB(file_path + '.sqlite3')
        else:
            gff_in = gffutils.create_db(file_path, dbfn=file_path + '.sqlite3')
        for feature in gff_in.featuretypes():
            logger.info("The number of %s in the gff file: %d", feature, gff_in.count_features_of_type(feature))
        return gff_in
    except Exception as e:
        logger.error("Failed to read the gff file: %s", file_path)
        raise e
