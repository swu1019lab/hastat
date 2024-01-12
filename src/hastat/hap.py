# -*- coding: utf-8 -*-
# @Time    : 2024/1/12 17:30
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : hap.py

import pandas as pd


class HapData(object):
    """
    A class for haplotype analysis
    """
    def __init__(self, geno_dataframe) -> None:
        """

        :param geno_dataframe: a dataframe contained genotype data
        """
        nrow, ncol = geno_dataframe.shape
        print("The number of output samples and loci:", ncol, nrow)
        # Seek haplotype of gene
        # access haplotypes and sequence
        hap_grouped = geno_dataframe.T.assign(size=1).groupby(geno_dataframe.index.to_list(), as_index=False)
        # haplotypes table, including name, genotype code and size
        # can be exported as fasta format
        hap_table = hap_grouped.sum()
        hap_table = hap_table.rename(index=lambda x: "Hap" + str(x + 1))
        hap_table.index.name = 'haplotypes'
        # haplotypes group
        # add haplotype labels, e.g. Hap1, Hap2, ..., HapN
        hap_ngroup = hap_grouped.ngroup()
        hap_ngroup.name = 'haplotypes'
        hap_ngroup = hap_ngroup.apply(lambda x: "Hap" + str(x + 1))
        # code number should transform into code character
        # i: chr, pos, ref, alt
        for i, s in hap_table.iloc[:, :-1].items():
            hap_table.loc[:, i] = s.map({0: i[2] + i[2], 1: i[2] + i[3], 2: i[3] + i[3]})
        hap_table.columns = hap_table.columns.droplevel(['ref', 'alt'])
        # access haplotypes group, number, counts, table and dataframe
        self.hap_ngroup = hap_ngroup.reset_index()
        self.hap_num = len(hap_grouped)
        self.hap_count = hap_ngroup.value_counts()
        self.hap_table = hap_table
        self.hap_dataframe = geno_dataframe.T.assign(haplotypes=hap_ngroup)

    def count_hap_of_groups(self, sample_groups=None):
        """
        Count haplotypes distribution within different samples source
        A sample group Dataframe with two columns only: samples and groups
        | samples | groups |
        | ------- | ------ |
        | ------- | ------ |

        :param sample_groups: add additional sample groups for haplotypes distribution count
        """
        if sample_groups is not None:
            hap_table0 = self.hap_ngroup.merge(sample_groups, how='inner', on='samples').groupby(
                ['haplotypes', 'groups']).count().unstack(fill_value=0)
            return hap_table0.loc[:, 'samples'].T
        else:
            print("Not found sample groups!")

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
            print("Getting the samples belong to {} from haplotypes groups".format(hap))
            return self.hap_ngroup.query('haplotypes == @hap').samples.to_list()

    def get_hap_nlargest(self, n=3, keep='first'):
        """
        Get the largest N haplotypes

        :param n: top N haplotypes for fetching, default is 3
        :param keep: default is first
        :return: a list containing haplotypes name
        """
        return self.hap_count.nlargest(n=n, keep=keep).index.tolist()

    def get_hap_ngroup(self):
        """
        Get all haplotypes and count

        :return: a dataframe with two columns, which can be used as p-value calculation
        """
        return self.hap_ngroup

    def get_hap_table(self, n=-1, sample_groups=None):
        """
        Get a haplotypes table

        :return: a dataframe containing detailed haplotypes information, which can be converted to other format
        """
        if sample_groups is not None:
            hap_table0 = self.hap_ngroup.merge(sample_groups, how='inner', on='samples').groupby(
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
