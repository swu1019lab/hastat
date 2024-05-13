# -*- coding: utf-8 -*-
# @Time    : 2024/1/12 17:30
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : hap.py

import pandas as pd
import numpy as np
from hastat.hastat import logger


class HapData(object):
    """
    A class for haplotype analysis
    """

    def __init__(self, df: pd.DataFrame) -> None:
        """
        Initialize the HapData object

        :param df: a dataframe contained genotype data. The format should be:
            1. the columns are samples;
            2. the rows are multi-indexed SNPs with chr, pos, ref and alt;
            3. the values are genotypes with 0, 1 and 2
        """
        # logger.info("Initializing the HapData object")
        # # access gene haplotype data
        # hap_grouped = df.T.assign(size=1).groupby(df.index.to_list(), as_index=False)
        #
        # # get haplotypes table, including name, genotype code and size
        # hap_table = hap_grouped.sum().sort_values(by='size', ascending=False).reset_index()
        # hap_table = hap_table.rename(index=lambda x: "Hap" + str(x + 1))
        # hap_table.index.name = 'haplotypes'
        #
        # # get haplotypes group for each sample
        # hap_dict = hap_table.reset_index().set_index('index')['haplotypes'].to_dict()
        # hap_groups = hap_grouped.ngroup()
        # hap_groups.name = 'haplotypes'
        # hap_groups = hap_groups.map(hap_dict)
        #
        # # transform code number into code character
        # # i: chr, pos, ref, alt
        # hap_table = hap_table.drop(columns='index')
        # print(hap_table.columns)
        # for i, s in hap_table.iloc[:, :-1].items():
        #     hap_table.loc[:, i] = s.map({0: i[2] + i[2], 1: i[2] + i[3], 2: i[3] + i[3]})
        # hap_table.columns = hap_table.columns.droplevel(['ref', 'alt'])
        # print(hap_table.columns)
        #
        # # access haplotypes group, number, counts, table and dataframe
        # self.hap_groups = hap_groups.reset_index()
        # self.hap_num = hap_grouped.ngroups
        # self.hap_count = hap_groups.value_counts()
        # self.hap_table = hap_table
        # self.hap_dataframe = df.T.assign(haplotypes=hap_groups)
        logger.info("Initializing the HapData object")

        # get unique haplotypes and their counts
        hap, counts = np.unique(df.T.values, axis=0, return_counts=True)
        _, indices = np.unique(df.T.values, axis=0, return_inverse=True)

        # Access haplotype data
        # hap_grouped = df.T.assign(size=1).groupby(df.index.to_list(), as_index=False)

        # haplotypes table, including name, genotype code and size
        hap_table = pd.DataFrame(hap, columns=df.index)
        hap_table['size'] = counts
        hap_table.index = "Hap" + (hap_table.index + 1).astype(str)
        hap_table.index.name = 'haplotypes'

        # haplotypes group for each sample
        # add haplotype labels, e.g. Hap1, Hap2, ..., HapN
        hap_groups = pd.Series(indices, index=df.columns)
        hap_groups = "Hap" + (hap_groups + 1).astype(str)
        hap_groups.name = 'haplotypes'

        # code number should transform into code character
        # i: chr, pos, ref, alt
        for i, s in hap_table.iloc[:, :-1].items():
            ref, alt = i[2], i[3]
            hap_table.loc[:, i] = s.map({0: ref + ref, 1: ref + alt, 2: alt + alt})
        hap_table.columns = hap_table.columns.droplevel(['ref', 'alt'])

        # access haplotypes group, number, counts, table and dataframe
        self.hap_groups = hap_groups.reset_index()
        self.hap_num = len(hap)
        self.hap_count = hap_groups.value_counts()
        self.hap_table = hap_table
        self.hap_dataframe = df.T.assign(haplotypes=hap_groups)
        # Convert haplotypes to a categorical data type to save memory and improve performance
        self.hap_dataframe['haplotypes'] = self.hap_dataframe['haplotypes'].astype('category')

    def get_hap_freq(self, sample_groups=None):
        """
        Count haplotypes distribution within different samples source
        A sample group Dataframe with two columns only: samples and groups
        | samples | groups |
        | ------- | ------ |
        | ------- | ------ |

        :param sample_groups: add additional sample groups for haplotypes distribution count
        """
        if sample_groups is not None:
            hap_table0 = self.hap_groups.merge(sample_groups, how='inner', on='samples').groupby(
                ['haplotypes', 'groups']).count().unstack(fill_value=0)
            return hap_table0.loc[:, 'samples'].T
        else:
            return self.hap_count

    def get_snps_data(self):
        """
        Get SNPs data

        :return: a dataframe containing snps data
        """
        return self.hap_table.columns[:-1]

    def get_samples_of_hap(self, hap=None):
        """
        Get samples belong to target haplotype

        :param hap: haplotype name
        :return: a list containing samples from target haplotype
        """
        if hap is not None:
            logger.info("Getting the samples belong to {} from haplotypes groups".format(hap))
            return self.hap_groups.query('haplotypes == @hap').samples.to_list()

    def get_hap_nlargest(self, n=3, keep='first'):
        """
        Get the largest N haplotypes

        :param n: top N haplotypes for fetching, default is 3
        :param keep: default is first
        :return: a list containing haplotypes name
        """
        return self.hap_count.nlargest(n=n, keep=keep).index.tolist()

    def get_hap_groups(self):
        """
        Get all haplotypes and count

        :return: a dataframe with two columns, which can be used as p-value calculation
        """
        return self.hap_groups

    def get_hap_table(self, n=-1, sample_groups=None):
        """
        Get a haplotypes table

        :return: a dataframe containing detailed haplotypes information, which can be converted to other format
        """
        if sample_groups is not None:
            hap_table0 = self.hap_groups.merge(sample_groups, how='inner', on='samples').groupby(
                ['haplotypes', 'groups']).count().unstack(fill_value=0)
            self.hap_table = self.hap_table.merge(hap_table0, left_index=True, right_index=True)
        if n > 0:
            # Return the first n rows ordered by columns in descending order.
            return self.hap_table.nlargest(n, 'size')
        return self.hap_table

    def get_hap_dataframe(self):
        """
         Get a haplotypes dataframe

        :return: a dataframe contained genotype and haplotype information
        """
        return self.hap_dataframe

    def to_fasta(self):
        """
        Export haplotypes data into fasta format

        :return: None
        """
        for hap, seq in self.hap_table.iterrows():
            print(f">{hap}\n{seq[:-1].str.cat()}")

    def to_phylip(self):
        """
        Export haplotypes data into phylip format

        :return: None
        """
        nrow, ncol = self.hap_table.shape
        print(f" {nrow} {ncol - 1}")
        for hap, seq in self.hap_table.iterrows():
            print("{:<10s}{}".format(hap, seq[:-1].str.cat()))

    def to_excel(self, file_path='haplotypes.data.xlsx'):
        """
        Export haplotypes data into Excel file

        :param file_path: the path to save haplotypes data, default is haplotypes.data.xlsx
        :return: None
        """
        with pd.ExcelWriter(file_path) as writer:
            self.hap_count.to_excel(writer, sheet_name='count')
            self.hap_table.to_excel(writer, sheet_name='table')
            self.hap_dataframe.to_excel(writer, sheet_name='dataframe')
