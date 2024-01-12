# -*- coding: utf-8 -*-
# @Time    : 2024/1/12 17:30
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : stat.py

from collections import defaultdict
import pandas as pd
from statsmodels.stats.oneway import anova_oneway
from statsmodels.sandbox.stats.multicomp import MultiComparison


class StatData(object):
    """
    a class for oneway anova and multi-comparison
    """

    def __init__(self, pheno, groups):
        """_summary_

        Args:
            pheno (DataFrame): a two-columns at least DataFrame contained "samples" and "other traits"
            groups (DataFrame): a two-columns DataFrame contained "samples" and "haplotypes"
        """
        self.phap_data = pheno.merge(groups, how='inner', on='samples')
        self.phap_data.set_index(['samples', 'haplotypes'], inplace=True)

        self._stat_des = self.phap_data.groupby(level='haplotypes').describe().fillna(0)
        self._stat_anova = defaultdict(list)

        # Oneway anova and multiple comparisons results
        for phe in self.phap_data.columns:
            # filter haplotype with na value and std equals to zero
            # (which means no difference between means of all groups)
            pass_hap = self._stat_des.loc[:, (phe,)].query('std>0').index.to_list()

            # 2 or more groups required for multiple comparisons
            if len(pass_hap) > 1:
                pass_phap_data = self.phap_data.loc[(slice(None), pass_hap), phe].dropna().reset_index()
                # oneway anova
                anova_res = anova_oneway(pass_phap_data.iloc[:, 2], groups=pass_phap_data.iloc[:, 1])
                # print(phe, anova_res.statistic, anova_res.pvalue)
                self._stat_anova['pvalue'].append(anova_res.pvalue)
                # Tests for multiple comparisons
                multi_comp = MultiComparison(pass_phap_data.iloc[:, 2], groups=pass_phap_data.iloc[:, 1])
                multi_comp_tbl = multi_comp.tukeyhsd(alpha=0.05).summary()
                midx = pd.MultiIndex.from_product([[phe], multi_comp_tbl.data[0]])
                multi_comp_df = pd.DataFrame(data=multi_comp_tbl.data[1:], columns=midx)
                self._stat_anova['comparisons'].append(multi_comp_df)
            else:
                print("Warning: 2 or more groups required for oneway anova!")
                print("\tFailed to perform multiple comparisons for phenotype {}".format(phe))
                self._stat_anova['pvalue'].append("")

    def get_basic_data(self):
        # Generate descriptive statistics of each haplotype for each phenotype, like min, max, mean, std and so on.
        return self._stat_des

    def get_anova_data(self):
        # Oneway anova and multiple comparisons results
        return self._stat_anova

    def to_excel(self, file_path='haplotypes.stats.xlsx'):
        # Export results into Excel file
        with pd.ExcelWriter(file_path) as writer:
            self._stat_des.to_excel(writer, sheet_name='basic')
            if len(self._stat_anova['pvalue']) > 0:
                stat_multi_pvalue = pd.DataFrame([self._stat_anova['pvalue']], columns=self.phap_data.columns)
                stat_multi_pvalue.to_excel(writer, sheet_name='pvalue')
            if len(self._stat_anova['comparisons']) > 0:
                stat_multi_comp = pd.concat(self._stat_anova['comparisons'], axis=1)
                stat_multi_comp.to_excel(writer, sheet_name='comparisons')
