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
from hastat.viz._export import save_figure
from matplotlib.colors import LinearSegmentedColormap
from hastat.dataio.table import read_sample_csv, normalize_sample_table, select_haplotypes, merge_samples, grouped_counts


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
        self.add_hap(read_sample_csv(data['sample_hap']), self.config['plot'].get('haplotypes'))
        self.add_group(read_sample_csv(data['sample_group']))

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
        self.sample_hap = select_haplotypes(normalize_sample_table(df, 'haplotypes'), haplotypes)

    def add_group(self, df: pd.DataFrame):
        """
        Add a dataframe contained group data

        :param df: a dataframe contained group data, the format should be:
            1. two columns at least
            2. 1st column is sample names
            3. other columns are group names
        :return:
        """
        self.sample_group = normalize_sample_table(df, 'groups')

    def merge_data(self):
        if self.sample_hap is None or self.sample_group is None:
            raise ValueError("hap_data or group_data is None")
        self.data = merge_samples(self.sample_hap, self.sample_group)

    def get_plot_data(self, name: str, calc_percentage: bool = False):
        self.merge_data()
        return grouped_counts(self.data, name, calc_percentage)

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
            ax.text(i * 1, 0.5 + 0.48, col, ha='center', va='bottom', fontsize=12, rotation=0)
            ax.text(i * 1, 0.5 - 0.48, f"{data[col].sum()}", ha='center', va='top', fontsize=11)
            ax.axis('equal')
        ax.legend(data.index, bbox_to_anchor=(1, 0, 0.5, 1), loc='center left', frameon=False,
                  fontsize=11, labelspacing=.7)
        save_figure(ax.figure, self.config['plot']['save_fig'], dpi=300, bbox_inches='tight')
