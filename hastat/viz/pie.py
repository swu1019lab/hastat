# -*- coding: utf-8 -*-
# @Time    : 2024/9/30 18:42
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : pie.py


import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import axes


class HapPie(object):
    def __init__(self):
        """
        Initialize HapPie object
        """
        self.hap_group = None
        self.sample_group = None
        self.data = None

    def add_hap(self, df: pd.DataFrame, labels: list = None):
        """
        Add a dataframe contained haplotype data

        :param df: a dataframe contained haplotype data, the format should be:
            1. two columns
            2. 1st column is sample names
            3. 2nd column is haplotype groups
        :param labels: a list contained haplotype group names to be analyzed
        :return:
        """
        if not isinstance(df, pd.DataFrame):
            raise ValueError("df is not a DataFrame")
        df.columns = ['sample', 'haplotypes']
        if labels is None:
            self.hap_group = df
        else:
            self.hap_group = df[df['haplotypes'].isin(labels)]

    def add_group(self, df: pd.DataFrame):
        """
        Add a dataframe contained group data

        :param df: a dataframe contained group data, the format should be:
            1. two columns at least
            2. 1st column is sample names
            3. other columns are group names
        :return:
        """
        if not isinstance(df, pd.DataFrame):
            raise ValueError("df is not a DataFrame")
        df.columns = ['sample'] + df.columns[1:].tolist()
        self.sample_group = df

    def merge_data(self):
        if self.hap_group is None or self.sample_group is None:
            raise ValueError("hap_data or group_data is None")
        self.data = self.hap_group.merge(self.sample_group, how='inner', on='sample')

    def get_plot_data(self, name: str, calc_percentage: bool = False):
        self.merge_data()
        data = self.data.loc[:, ['sample', 'haplotypes', name]].groupby(['haplotypes', name]).size().unstack(
            fill_value=0)
        if calc_percentage:
            return data.div(data.sum())
        return data

    def plot(self, ax: axes.Axes, name: str = None, colors: list = None, **kwargs):
        """
        Plot pie chart for haplotype data

        :param ax: axes object
        :param name: group name to be analyzed
        :param colors: colors for each haplotype
        :param kwargs: other parameters for pie chart
        :return:
        """
        if name is None:
            name = self.sample_group.columns[1]
        data = self.get_plot_data(name, calc_percentage=False)
        if colors is None:
            colors = mpl.colormaps['Blues'](np.linspace(0.2, 0.7, len(data.index)))
        # plot pie chart for each column
        for i, col in enumerate(data.columns):
            ax.pie(data[col], labels=data.index, labeldistance=None,
                   colors=colors, center=(i * 1, 0.5), radius=0.4,
                   wedgeprops={"linewidth": 1, "edgecolor": "white"},
                   **kwargs)
            ax.text(i * 1, 0.5 + 0.5, col, ha='center', va='bottom', fontsize=8, rotation=0)
            ax.text(i * 1, 0.5 - 0.5, f"{data[col].sum()}", ha='center', va='top', fontsize=8)
            ax.axis('equal')
        ax.legend(data.index, bbox_to_anchor=(1, 0, 0.5, 1), loc='center left', frameon=False)
