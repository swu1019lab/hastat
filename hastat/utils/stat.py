# -*- coding: utf-8 -*-
# @Time    : 2024/1/12 17:30
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : stat.py

from collections import defaultdict
import pandas as pd
from statsmodels.stats.oneway import anova_oneway
from statsmodels.stats.weightstats import ttest_ind
from statsmodels.sandbox.stats.multicomp import MultiComparison
from hastat.hastat import logger


class HapAnovaTest(object):
    """
    a class for oneway anova and multi-comparison
    """
    def __init__(self, pheno: pd.DataFrame, groups: pd.DataFrame, min_hap_size=10, annotate='gene', method='TukeyHSD'):
        """ Initialize the class with phenotype and haplotype groups

        Args:
            pheno (DataFrame): a two-columns at least DataFrame contained "samples" and other traits name
            groups (DataFrame): a two-columns DataFrame contained "samples" and "haplotypes"
            min_hap_size (int): the minimum sample size for each haplotype, default is 10
            annotate (str): the annotation of haplotypes, default is "gene"
            method (str): the method for multiple comparisons, default is 'TukeyHSD', also can be AllPairTest
        """
        # Filter haplotypes with size less than min_hap_size (only keep major haplotypes)
        groups = groups[groups['haplotypes'].map(groups['haplotypes'].value_counts() >= min_hap_size)]
        if groups.empty or groups['haplotypes'].nunique() < 2:
            raise ValueError(
                "Two or more haplotypes required for one-way ANOVA! "
                "You can try to reduce the min_hap_size.")
        # Merge phenotype and haplotype groups
        self.phap_data = pheno.merge(groups, how='inner', on='samples')
        self.phap_data.set_index(['samples', 'haplotypes'], inplace=True)

        # Generate descriptive statistics of each haplotype for each phenotype, like min, max, mean, std and so on.
        self._stat_des = self.phap_data.groupby(level='haplotypes').describe().fillna(0)
        self._stat_anova = defaultdict(list)

        # Oneway anova and multiple comparisons results
        for phe in self.phap_data.columns:
            self.perform_anova_and_comparisons(phe, annotate=annotate, method=method)

    def perform_anova_and_comparisons(self, phe: str, annotate='gene', method='TukeyHSD'):
        """
        Perform oneway anova and multiple comparisons for each phenotype using statsmodels

        :param phe: phenotype name
        :param annotate: the annotation of haplotypes
        :param method: the method for multiple comparisons, default is 'TukeyHSD', also can be AllPairTest
        :return: None
        """
        # filter haplotype with na value and std equals to zero
        # (which means no difference between means of all groups)
        pass_hap = self._stat_des.loc[:, (phe,)].query('std>0').index.to_list()

        # 2 or more groups required for multiple comparisons
        if len(pass_hap) > 1:
            pass_phap_data = self.phap_data.loc[(slice(None), pass_hap), phe].dropna().reset_index()
            # oneway anova
            anova_res = anova_oneway(pass_phap_data.iloc[:, 2], groups=pass_phap_data.iloc[:, 1])

            # Tests for multiple comparisons
            multi_comp = MultiComparison(pass_phap_data.iloc[:, 2], groups=pass_phap_data.iloc[:, 1])
            if method == 'AllPairTest':
                multi_comp_tbl = multi_comp.allpairtest(ttest_ind, alpha=0.05, method="bonf")[0]
            else:
                multi_comp_tbl = multi_comp.tukeyhsd(alpha=0.05).summary()
            # index = pd.MultiIndex.from_product([[phe], multi_comp_tbl.data[0]])
            # Add the number of haplotypes and phenotype mean for each group, optional statistics like std, min, max
            pass_phap_count = self._stat_des.loc[:, (phe, 'count')].to_dict()
            pass_phap_mean = self._stat_des.loc[:, (phe, 'mean')].to_dict()
            multi_comp_df = pd.DataFrame(data=multi_comp_tbl.data[1:], columns=multi_comp_tbl.data[0])
            multi_comp_df = multi_comp_df.assign(pheno=phe,
                                                 annotate=annotate,
                                                 anova=anova_res.pvalue,
                                                 count1=lambda x: x['group1'].map(pass_phap_count),
                                                 count2=lambda x: x['group2'].map(pass_phap_count),
                                                 mean1=lambda x: x['group1'].map(pass_phap_mean),
                                                 mean2=lambda x: x['group2'].map(pass_phap_mean))

            self._stat_anova['comparisons'].append(multi_comp_df)
        else:
            logger.warning(
                f"Two or more groups required for one-way ANOVA! Failed to "
                f"perform multiple comparisons for phenotype {phe}")

    @property
    def basic_data(self):
        return self._stat_des

    @property
    def anova_data(self):
        return self._stat_anova
