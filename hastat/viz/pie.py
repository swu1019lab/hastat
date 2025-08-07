# -*- coding: utf-8 -*-
# @Time    : 2024/9/30 18:42
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : pie.py


import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import axes
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


class HapPie(object):
    def __init__(self, config: dict = None):
        """
        Initialize HapPie object
        """
        if config is None:
            raise ValueError("config is None")
        self.config = config
        self.sample_hap = None
        self.sample_group = None
        self.data = None

    def add_data(self, data: dict):
        self.add_hap(pd.read_csv(data['sample_hap']), self.config['plot']['haplotypes'])
        self.add_group(pd.read_csv(data['sample_group']))

    def add_hap(self, df: pd.DataFrame, haplotypes: list = None):
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
        if haplotypes is None or len(haplotypes) == 0:
            self.sample_hap = df
        else:
            self.sample_hap = df[df['haplotypes'].isin(haplotypes)]

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
        if self.sample_hap is None or self.sample_group is None:
            raise ValueError("hap_data or group_data is None")
        self.data = self.sample_hap.merge(self.sample_group, how='inner', on='sample')

    def get_plot_data(self, name: str, calc_percentage: bool = False):
        self.merge_data()
        data = self.data.loc[:, ['sample', 'haplotypes', name]].groupby(['haplotypes', name]).size().unstack(
            fill_value=0)
        if calc_percentage:
            return data.div(data.sum())
        return data

    def plot(self, ax: axes.Axes = None):
        """
        Plot pie chart for haplotype data

        :param ax: axes object
        :return:
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(self.config['plot']['width'], self.config['plot']['height']))

        self.add_data(self.config['data'])

        name = self.sample_group.columns[self.config['plot']['group_index']]
        data = self.get_plot_data(name, calc_percentage=self.config['plot']['calc_percentage'])
        pie_cmap = LinearSegmentedColormap.from_list('pie_cmap', ['#C5504B', '#114F8B', '#FCE988', '#90CAEE'], N=100)
        pie_colors = pie_cmap(np.linspace(0, 1, len(data.index)))
        # plot pie chart for each column
        for i, col in enumerate(data.columns):
            ax.pie(data[col], labels=data.index, labeldistance=None,
                   colors=pie_colors, center=(i * 1, 0.5), radius=0.4,
                   wedgeprops={"linewidth": 1, "edgecolor": "white"})
            ax.text(i * 1, 0.5 + 0.5, col, ha='center', va='bottom', fontsize=8, rotation=0)
            ax.text(i * 1, 0.5 - 0.5, f"{data[col].sum()}", ha='center', va='top', fontsize=8)
            ax.axis('equal')
        ax.legend(data.index, bbox_to_anchor=(1, 0, 0.5, 1), loc='center left', frameon=False)
        plt.savefig(self.config['plot']['save_fig'], dpi=300, bbox_inches='tight')
