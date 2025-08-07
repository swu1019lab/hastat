# -*- coding: utf-8 -*-
# @Time    : 2024/10/8 16:01
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : box.py

import pandas as pd
import numpy as np
from matplotlib import axes
from collections import Counter
import matplotlib.pyplot as plt
from scipy import stats
from hastat.log.logger import logger
from matplotlib.colors import LinearSegmentedColormap


class HapBox(object):
    def __init__(self, config: dict = None):
        if config is None:
            raise ValueError("config is None")
        self.config = config

        self.sample_hap = None
        self.sample_phe = None
        self.comparisons = None
        self.method = None
        self.step_size = None
        self.data = None

    def add_data(self, data: dict):
        self.add_hap(pd.read_csv(data['sample_hap']), self.config['plot']['haplotypes'])
        self.add_phe(pd.read_csv(data['sample_phe']))

    def add_hap(self, df: pd.DataFrame, haplotypes: list = None):
        """
        Add a dataframe contained haplotype data

        :param df: a dataframe contained haplotype data, the format should be:
            1. two columns
            2. 1st column is sample names
            3. 2nd column is haplotype groups
        :param haplotypes: a list contained haplotype group names to be analyzed
        :return:
        """
        if not isinstance(df, pd.DataFrame):
            raise ValueError("df is not a DataFrame")
        df.columns = ['sample', 'haplotypes']
        if haplotypes is None or len(haplotypes) == 0:
            self.sample_hap = df
        else:
            self.sample_hap = df[df['haplotypes'].isin(haplotypes)]

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
            self.sample_phe = df
        else:
            self.sample_phe = df.loc[:, ['sample'] + name]

    def add_comparisons(self, comparisons: list):
        """
        Add comparisons to perform significant test
        :param comparisons: a list contained comparisons, e.g. [('A', 'B'), ('A', 'C')].
        :return:
        """
        self.comparisons = comparisons

    def add_method(self, method: str = 't-test'):
        """
        Add statistical test method

        :param method: statistical test method, e.g. 't-test', 'mannwhitneyu', 'welch'.
        :return:
        """
        self.method = method

    def add_step_size(self, step_size: int or list):
        """
        Add step size to the plot

        :param step_size: a list contained step size for each comparison.
        :return:
        """
        if not isinstance(step_size, list):
            self.step_size = [step_size] * len(self.comparisons)
        else:
            self.step_size = step_size

    def perform_statistical_test(self, group1_data, group2_data):
        """
        Perform statistical test between two groups
        
        :param group1_data: data for group 1
        :param group2_data: data for group 2
        :return: p-value
        """
        if self.method == 't-test':
            # Perform two-sample t-test
            _, p_value = stats.ttest_ind(group1_data, group2_data)
        elif self.method == 'mannwhitneyu':
            # Perform Mann-Whitney U test
            _, p_value = stats.mannwhitneyu(group1_data, group2_data)
        elif self.method == 'welch':
            # Perform Welch's t-test (unequal variances)
            _, p_value = stats.ttest_ind(group1_data, group2_data, equal_var=False)
        else:
            raise ValueError(f"Unknown statistical method: {self.method}")
        
        return p_value

    def get_significance_symbol(self, p_value):
        """
        Get significance symbol based on p-value
        
        :param p_value: p-value from statistical test
        :return: significance symbol
        """
        if p_value < 0.001:
            return '***'
        elif p_value < 0.01:
            return '**'
        elif p_value < 0.05:
            return '*'
        else:
            return 'ns'

    def merge_data(self):
        if self.sample_hap is None or self.sample_phe is None:
            raise ValueError("hap_data or sample_phe is None")
        self.data = self.sample_hap.merge(self.sample_phe, how='inner', on='sample')

    def get_plot_data(self, name: str):
        self.merge_data()
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

    def plot(self, ax: axes.Axes = None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(self.config['plot']['width'], self.config['plot']['height']))

        self.add_data(self.config['data'])
        
        # Set comparisons and method if provided
        if self.config['plot']['comparisons']:
            self.add_comparisons(self.config['plot']['comparisons'])
            self.add_method(self.config['plot']['method'])
            if self.config['plot']['step_size']:
                self.add_step_size(self.config['plot']['step_size'])

        name = self.sample_phe.columns[self.config['plot']['phe_index']]
        xdata, ydata, data, labels = self.get_plot_data(name)

        bp = ax.boxplot(list(data.values()),
                        labels=labels,
                        widths=0.5,
                        showfliers=False,
                        # fill patch with color
                        patch_artist=True
                        )
        box_cmap = LinearSegmentedColormap.from_list('box_cmap', ['#C5504B', '#114F8B', '#FCE988', '#90CAEE'], N=100)
        box_colors = box_cmap(np.linspace(0, 1, len(bp['boxes'])))
        for p, c in zip(bp['boxes'], box_colors):
            p.set_facecolor(c)
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

        # add comparisons with automatic statistical testing
        group_count = Counter()
        if self.comparisons is not None:
            for i, (group1, group2) in enumerate(self.comparisons):
                if group1 in data and group2 in data:
                    y1 = np.array(data[group1])
                    y2 = np.array(data[group2])
                    x1 = xdata.categories.get_loc(group1) + 1
                    x2 = xdata.categories.get_loc(group2) + 1
                    group_count[group1] += 1
                    group_count[group2] += 1
                    # calculate mean and std of two groups
                    logger.info(f'group {group1}: {np.round(np.mean(y1), 2)}±{np.round(np.std(y1), 2)} vs group {group2}: {np.round(np.mean(y2), 2)}±{np.round(np.std(y2), 2)}')

                    # Perform statistical test
                    p_value = self.perform_statistical_test(y1, y2)
                    sig_symbol = self.get_significance_symbol(p_value)
                    
                    # add significant line or bracket
                    max_n = max(group_count[group1], group_count[group2])
                    max_y = np.max(np.concatenate((y1, y2))) * (1 + 0.02 * max_n)
                    if self.step_size is not None and i < len(self.step_size):
                        max_y += self.step_size[i]
                    ax.plot([x1, x2], [max_y, max_y], marker=3, markersize=5, lw=1, c='k')
                    # add significant symbol
                    ax.text((x1 + x2) / 2, max_y, sig_symbol, ha='center', va='bottom')

        # set the axes ranges
        ax.set(
            axisbelow=True,
            xlim=(0.5, len(bp.get('boxes')) + 0.5),
            ylim=(None, None),
        )
        ax.set_ylabel(name)
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.spines[["top", "right"]].set_visible(False)
        ax.yaxis.grid(
            True,
            linestyle='--',
            which='major',
            color='lightgrey',
            alpha=.5
        )
        plt.savefig(self.config['plot']['save_fig'], dpi=300, bbox_inches='tight')
