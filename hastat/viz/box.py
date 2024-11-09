# -*- coding: utf-8 -*-
# @Time    : 2024/10/8 16:01
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : box.py

import pandas as pd
import numpy as np
from itertools import cycle
from matplotlib import axes
from collections import Counter


class HapBox(object):
    def __init__(self):
        self.hap_group = None
        self.phe_data = None
        self.comparisons = None
        self.sig_symbol = None
        self.step_size = None
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

    def add_phe(self, df: pd.DataFrame, name: str or list = None):
        """
        Add a dataframe contained phenotype data

        :param df: a dataframe contained phenotype data, the format should be:
            1. two columns at least
            2. 1st column is sample names
            3. other columns are phenotype values
        :param name: the name of phenotype to be added
        :return:
        """
        if not isinstance(df, pd.DataFrame):
            raise ValueError("df is not a DataFrame")

        df.columns = ['sample'] + df.columns[1:].tolist()
        if name is None:
            self.phe_data = df
        else:
            self.phe_data = df.loc[:, ['sample'] + name]

    def add_comparisons(self, comparisons: list):
        """
        Add comparisons to perform significant test
        :param comparisons: a list contained comparisons, e.g. [('A', 'B'), ('A', 'C')].
        :return:
        """
        self.comparisons = comparisons

    def add_sig_symbol(self, sig_symbol: list = None):
        """
        Add significant symbol to the plot

        :param sig_symbol: a list contained significant symbol, e.g. ['*', '**', '***'].
        :return:
        """
        self.sig_symbol = sig_symbol

    def add_step_size(self, step_size: int or list):
        """
        Add step size to the plot

        :param step_size: a list contained step size for each comparison.
        :return:
        """
        if isinstance(step_size, int):
            self.step_size = [step_size] * len(self.comparisons)
        else:
            self.step_size = step_size

    def merge_data(self):
        if self.hap_group is None or self.phe_data is None:
            raise ValueError("hap_data or phe_data is None")
        self.data = self.hap_group.merge(self.phe_data, how='inner', on='sample')

    def get_plot_data(self, name: str):
        x0data = self.data['haplotypes'].to_numpy()
        y0data = self.data[name].to_numpy(dtype=float)

        # Process infinite or NaN data
        indices = np.nonzero(np.isfinite(y0data))[0]

        xdata = pd.Categorical(x0data[indices])
        ydata = y0data[indices]
        labels = xdata.categories

        # Add value according to the order of categories
        data = {k: [] for k in labels}
        for x, y in zip(xdata, ydata):
            data[x].append(y)
        return xdata, ydata, data, labels

    def plot(self, ax: axes.Axes, name: str):
        if ax is None:
            raise ValueError("ax is None")
        if name is None:
            raise ValueError("name is None")

        self.merge_data()
        xdata, ydata, data, labels = self.get_plot_data(name)
        bp = ax.boxplot(list(data.values()),
                        labels=labels,
                        widths=0.5,
                        showfliers=False,
                        # fill patch with color
                        patch_artist=True
                        )
        box_colors = cycle(['#498DCB', '#F9BEBF', '#747474', '#EE3424'])
        for p in ax.patches:
            p.set_facecolor(next(box_colors))
            p.set_edgecolor('k')

        # add points to boxplot
        offsets = np.zeros((len(ydata), 2))
        np.random.seed(0)
        jitter = np.random.uniform(low=-0.25, high=0.25, size=len(ydata))
        offsets[:, 0] = xdata.codes + 1 + jitter
        offsets[:, 1] = ydata

        pc = ax.scatter(xdata.codes + 1, ydata, s=1)
        pc.set_facecolor('gray')
        # set jitter for scatter
        pc.set_offsets(offsets)
        # move the scatters on top of the line.
        pc.set_zorder(2)

        # add comparisons
        group_count = Counter()
        if self.comparisons is not None and self.sig_symbol is not None:
            for i, (group1, group2) in enumerate(self.comparisons):
                y1 = np.array(data[group1])
                y2 = np.array(data[group2])
                x1 = xdata.categories.get_loc(group1) + 1
                x2 = xdata.categories.get_loc(group2) + 1
                group_count[group1] += 1
                group_count[group2] += 1
                # add significant symbol
                sig = self.sig_symbol[i]
                # add significant line or bracket
                max_n = max(group_count[group1], group_count[group2])
                max_y = np.max(np.concatenate((y1, y2))) * (1 + 0.02 * max_n)
                if self.step_size is not None:
                    max_y += self.step_size[i]
                ax.plot([x1, x2], [max_y, max_y], marker=3, markersize=5, lw=1, c='k')
                # add significant symbol
                ax.text((x1 + x2) / 2, max_y, sig, ha='center', va='bottom')

        # set the axes ranges
        ax.set(
            axisbelow=True,
            title=name,
            xlim=(0.5, len(bp.get('boxes')) + 0.5),
            ylim=(None, None),
        )
        ax.set_ylabel('Value')
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.spines[["top", "right"]].set_visible(False)
        ax.yaxis.grid(
            True,
            linestyle='--',
            which='major',
            color='lightgrey',
            alpha=.5
        )
