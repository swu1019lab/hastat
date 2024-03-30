# -*- coding: utf-8 -*-
# @Time    : 2024/1/12 17:30
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : geno.py

import pandas as pd
import numpy as np
import allel
from scipy.spatial.distance import squareform
from .hap import HapData


class GenoData(object):
    """
    A class for processing genotype data
    """

    def __init__(self, allele_data_set) -> None:
        """
        :param allele_data_set: a dict contained genotype data
        """
        self.allele_data_set = allele_data_set

    def get_geno_data(self, code='code1'):
        """
        Get genotype data

        :param code: code1 or code2, default code1, code1 is presented as 0/1/2,
        code2 is presented as [0, 0], [0, 1] and [1, 1]
        :return: a DataFrame contained genotype data
        """
        return pd.DataFrame.from_dict(self.allele_data_set[code], orient='tight')

    def get_hap_data(self):
        """
        Get haplotype data

        :return: a HapData object
        """
        return HapData(self.get_geno_data())

    def get_snp_pos(self):
        """
        Get SNP position

        :return: a numpy array contained SNP position
        """
        return self.get_geno_data().index.get_level_values('pos').to_numpy()

    def get_pi_data(self, size=1, step=1, sample_group=None):
        """
        Get pi data of target samples using sliding window method

        :param size: window size
        :param step: step size
        :param sample_group: a list contained sample names
        :return: a numpy array contained pi data
        """
        geno_data = self.get_geno_data()
        if sample_group is not None:
            data = geno_data.loc[:, sample_group].to_numpy()
        else:
            data = geno_data.to_numpy()
        h = allel.HaplotypeArray(data)
        ac = h.count_alleles()
        pos = geno_data.index.get_level_values('pos').to_list()
        pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac, size, step=step)
        return np.column_stack([windows, n_bases, counts, pi])

    def get_fst_data(self, sample_group1: list, sample_group2: list, size=1, step=1):
        """
        Get Fst data of target samples using sliding window method

        :param sample_group1: a list contained sample names
        :param sample_group2: a list contained sample names
        :param size: sliding window size
        :param step: step size
        :return: a numpy array contained Fst data
        """
        geno_data = self.get_geno_data('code2')
        # Returns -1 for unmatched values
        idx = geno_data.columns
        index_group1 = idx.get_indexer(idx.where(idx.isin(sample_group1)).dropna())
        index_group2 = idx.get_indexer(idx.where(idx.isin(sample_group2)).dropna())
        sub_pops = [index_group1, index_group2]
        pos = geno_data.index.get_level_values('pos').to_list()
        g = allel.GenotypeArray(geno_data.to_numpy().tolist())
        _ = np.seterr(divide='ignore')
        fst, windows, counts = allel.windowed_weir_cockerham_fst(pos, g, sub_pops, size, step=step)
        return np.column_stack([windows, counts, fst])

    def get_ld_data(self):
        """
        Get LD data of all samples

        :return: a numpy array contained LD data
        """
        geno_data = self.get_geno_data('code1')
        r = allel.rogers_huff_r(geno_data.to_numpy())
        return squareform(np.power(r, 2))

    def get_heterozygosity_data(self):
        """
        Get heterozygosity data of all samples

        :return: a numpy array contained heterozygosity data
        """
        geno_data = self.get_geno_data('code2')
        g = allel.GenotypeArray(geno_data.to_numpy().tolist())
        # the rate of observed heterozygosity
        het_o = allel.heterozygosity_observed(g)
        # the expected rate of heterozygosity
        af = g.count_alleles().to_frequencies()
        het_e = allel.heterozygosity_expected(af, ploidy=2)
        return np.column_stack([het_o, het_e])

    def get_inbreeding_coefficient(self):
        """
        Get inbreeding coefficient of all samples

        :return: a numpy array contained inbreeding coefficient
        """
        geno_data = self.get_geno_data('code2')
        g = allel.GenotypeArray(geno_data.to_numpy().tolist())
        return allel.inbreeding_coefficient(g)
