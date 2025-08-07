# -*- coding: utf-8 -*-
# @Time    : 2024/10/8 16:01
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : bar.py

import pandas as pd
import numpy as np
from matplotlib import axes
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


class HapBar(object):
    def __init__(self, config: dict = None):
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
        if ax is None:
            fig, ax = plt.subplots(figsize=(self.config['plot']['width'], self.config['plot']['height']))

        self.add_data(self.config['data'])

        name = self.sample_group.columns[self.config['plot']['group_index']]
        data = self.get_plot_data(name, calc_percentage=self.config['plot']['calc_percentage'])
        bottom = np.zeros(data.shape[1])
        x = data.columns.get_level_values(name).to_numpy()
        bar_cmap = LinearSegmentedColormap.from_list('bar_cmap', ['#C5504B', '#114F8B', '#FCE988', '#90CAEE'], N=100)
        bar_colors = bar_cmap(np.linspace(0, 1, len(data.index)))
        for i, (hap, row) in enumerate(data.iterrows()):
            ax.bar(x, row.to_numpy(), width=0.5, bottom=bottom, label=hap, color=bar_colors[i])
            bottom += row.to_numpy()

        # set the axes
        ax.set_xlabel(self.config['plot']['x_label'])
        ax.set_ylabel(self.config['plot']['y_label'])
        ax.spines[["top", "right"]].set_visible(False)
        ax.yaxis.grid(
            True,
            linestyle='--',
            which='major',
            color='lightgrey',
            alpha=.5
        )
        # Add legend
        ax.legend(loc='lower left', frameon=False,  ncol=len(data.index), bbox_to_anchor=(0., 1.02, 1., .102))

        plt.savefig(self.config['plot']['save_fig'], dpi=300, bbox_inches='tight')
