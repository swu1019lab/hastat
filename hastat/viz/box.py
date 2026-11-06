# -*- coding: utf-8 -*-
# @Time    : 2024/10/8 16:01
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : box.py

import pandas as pd
import numpy as np
from matplotlib import axes
from collections import Counter
from itertools import combinations
import matplotlib.pyplot as plt
from hastat.viz._export import plot_fonts, save_figure
from scipy import stats
from statsmodels.stats.oneway import anova_oneway
from hastat.log.logger import logger
from matplotlib.colors import LinearSegmentedColormap
from hastat.dataio.table import read_sample_csv, normalize_sample_table, select_haplotypes, merge_samples, parse_numeric_values
from hastat.utils.association import HaplotypeAssociation

BOX_METHODS = ('t-test', 'mannwhitneyu', 'welch', 'TukeyHSD', 'WelchBonferroni')


class HapBox(object):
    def __init__(self, config: dict = None):
        if config is None:
            raise ValueError("config is None")
        self.config = config

        self.sample_hap = None
        self.sample_phe = None
        self.comparisons = None
        self.method = 't-test'
        self.step_size = None
        self.data = None
        self._all_sample_hap = None
        self._test_data = None
        self.min_hap_size = config.get('plot', {}).get('min_hap_size', 2)
        if (not isinstance(self.min_hap_size, (int, np.integer)) or isinstance(self.min_hap_size, bool)
                or self.min_hap_size < 2):
            raise ValueError('min_hap_size must be an integer >= 2')

    def add_data(self, data: dict):
        self.add_hap(read_sample_csv(data['sample_hap']), self.config['plot'].get('haplotypes'))
        self.add_phe(read_sample_csv(data['sample_phe']))

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
        self._all_sample_hap = normalize_sample_table(df, 'haplotypes')
        self.sample_hap = select_haplotypes(self._all_sample_hap, haplotypes)

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
        df = normalize_sample_table(df, 'phenotypes')
        if name is None:
            self.sample_phe = df
        else:
            names = [name] if isinstance(name, str) else list(name)
            self.sample_phe = df.loc[:, ['sample'] + names].copy()

    def add_comparisons(self, comparisons: list):
        """
        Add comparisons to perform significant test
        :param comparisons: a list contained comparisons, e.g. [('A', 'B'), ('A', 'C')].
        :return:
        """
        pairs, seen = [], set()
        for pair in comparisons:
            if isinstance(pair, str) or len(pair) != 2:
                raise ValueError('Each comparison requires two haplotype names')
            a, b = (str(x).strip() for x in pair)
            key = tuple(sorted((a, b)))
            if not a or not b or a == b or key in seen:
                raise ValueError('Comparisons must be unique pairs of different haplotypes')
            pairs.append((a, b))
            seen.add(key)
        self.comparisons = pairs

    def add_method(self, method: str = 't-test'):
        """
        Add statistical test method

        :param method: statistical test method, e.g. 't-test', 'mannwhitneyu', 'welch'.
        :return:
        """
        if method not in BOX_METHODS:
            raise ValueError(f'Unknown statistical method: {method}')
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
        if self.method not in ('t-test', 'mannwhitneyu', 'welch'):
            raise ValueError('Use comparison_pvalues for family-based methods')
        group1_data = np.asarray(group1_data, dtype=float)
        group2_data = np.asarray(group2_data, dtype=float)
        group1_data = group1_data[np.isfinite(group1_data)]
        group2_data = group2_data[np.isfinite(group2_data)]
        if min(len(group1_data), len(group2_data)) < self.min_hap_size:
            return np.nan
        if self.method != 'mannwhitneyu' and group1_data.var() == 0 and group2_data.var() == 0:
            return np.nan
        if self.method == 't-test':
            # Perform two-sample t-test
            _, p_value = stats.ttest_ind(group1_data, group2_data, equal_var=True, alternative='two-sided')
        elif self.method == 'mannwhitneyu':
            # Perform Mann-Whitney U test
            _, p_value = stats.mannwhitneyu(group1_data, group2_data, alternative='two-sided', method='auto')
        elif self.method == 'welch':
            # Perform Welch's t-test (unequal variances)
            _, p_value = stats.ttest_ind(group1_data, group2_data, equal_var=False, alternative='two-sided')
        else:
            raise ValueError(f"Unknown statistical method: {self.method}")
        
        return p_value

    def comparison_pvalues(self, name):
        """Adjusted P values for the full input family, independent of displayed brackets."""
        self.get_plot_data(name)
        method = 'WelchBonferroni' if self.method == 'welch' else self.method
        if method in ('TukeyHSD', 'WelchBonferroni'):
            frame = self._test_data.rename(columns={'sample': 'samples'})
            result = HaplotypeAssociation(frame[['samples', name]], frame[['samples', 'haplotypes']],
                                  min_hap_size=self.min_hap_size, method=method)
            return {tuple(sorted((row['group1'], row['group2']))): row['p-adj']
                    for row in result.comparisons.to_dict('records')}
        arrays = {label: part[name].dropna().to_numpy(dtype=float)
                  for label, part in self._test_data.groupby('haplotypes', observed=True)}
        arrays = {label: values for label, values in arrays.items() if len(values) >= self.min_hap_size}
        pairs = list(combinations(sorted(arrays), 2))
        adjusted = {}
        for a, b in pairs:
            p = self.perform_statistical_test(arrays[a], arrays[b])
            adjusted[a, b] = min(1.0, len(pairs) * p) if np.isfinite(p) else np.nan
        return adjusted

    def omnibus_test(self, name):
        """Test the same eligible input family used by pairwise comparisons.

        Parametric methods use the existing ordinary/Welch ANOVA convention.
        Mann-Whitney comparisons use the tie-corrected Kruskal-Wallis H test
        (Kruskal and Wallis, 1952, doi:10.1080/01621459.1952.10483441), with
        its chi-square approximation restricted to at least five values/group.
        """
        self.get_plot_data(name)
        arrays = [part[name].dropna().to_numpy(dtype=float)
                  for _, part in self._test_data.groupby('haplotypes', observed=True)]
        arrays = [a for a in arrays if len(a) >= self.min_hap_size]
        if len(arrays) < 2:
            return None
        if self.method == 'mannwhitneyu':
            if min(map(len, arrays)) < 5:
                logger.warning('Kruskal-Wallis annotation requires at least five observations per eligible group')
                return None
            if np.ptp(np.concatenate(arrays)) == 0:
                return None
            result = stats.kruskal(*arrays)
            label = 'Kruskal-Wallis'
        else:
            unequal = self.method in ('welch', 'WelchBonferroni')
            variances = np.array([a.var(ddof=1) for a in arrays])
            if (unequal and np.any(variances == 0)) or np.all(variances == 0):
                return None
            result = anova_oneway(arrays, use_var='unequal' if unequal else 'equal',
                                  welch_correction=True)
            label = 'Welch ANOVA' if unequal else 'ANOVA'
        if not np.isfinite(result.pvalue):
            return None
        return dict(test=label, pvalue=float(result.pvalue), groups=len(arrays),
                    n=sum(map(len, arrays)))

    @staticmethod
    def _p_label(pvalue):
        if pvalue == 0:
            return r'$P < 2.23\times10^{-308}$'
        if pvalue < .001:
            mantissa, exponent = f'{pvalue:.2e}'.split('e')
            return rf'$P = {mantissa}\times10^{{{int(exponent)}}}$'
        return rf'$P = {pvalue:.3g}$'

    def get_significance_symbol(self, p_value):
        """
        Get significance symbol based on p-value
        
        :param p_value: p-value from statistical test
        :return: significance symbol
        """
        if p_value is None or not np.isfinite(p_value) or not 0 <= p_value <= 1:
            return None
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
        merged = merge_samples(self._all_sample_hap, self.sample_phe)
        self._test_data = merged
        self.data = merged.loc[merged['sample'].isin(self.sample_hap['sample'])].copy()

    def get_plot_data(self, name: str):
        self.merge_data()
        if name in ('sample', 'haplotypes') or name not in self.sample_phe:
            raise ValueError('Select a phenotype column, not the sample column')
        self._test_data[name] = parse_numeric_values(self._test_data[name])
        self.data[name] = parse_numeric_values(self.data[name])
        x0data = self.data['haplotypes'].to_numpy()
        y0data = self.data[name].to_numpy(dtype=float)

        # Process infinite or NaN data
        indices = np.nonzero(np.isfinite(y0data))[0]
        if not len(indices):
            raise ValueError('No finite phenotype observations to plot')

        xdata = pd.Categorical(x0data[indices])
        ydata = y0data[indices]
        labels = xdata.categories

        # Add value according to the order of categories
        data = {k: [] for k in labels}
        for x, y in zip(xdata, ydata):
            data[x].append(y)
        return xdata, ydata, data, labels

    @plot_fonts()
    def plot(self, ax: axes.Axes = None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(self.config['plot']['width'], self.config['plot']['height']))

        self.add_data(self.config['data'])
        
        self.add_method(self.config['plot'].get('method', 't-test'))
        # Select comparison brackets independently of the overall test.
        self.comparisons = None
        if self.config['plot'].get('comparisons'):
            self.add_comparisons(self.config['plot']['comparisons'])
            if self.config['plot'].get('step_size'):
                self.add_step_size(self.config['plot']['step_size'])

        name = self.sample_phe.columns[self.config['plot']['phe_index']]
        xdata, ydata, data, labels = self.get_plot_data(name)
        rotation = self.config['plot'].get('x_rotation', 0)
        if not isinstance(rotation, (int, float, np.number)) or isinstance(rotation, bool) or not np.isfinite(rotation):
            raise ValueError('x_rotation must be a finite angle in degrees')
        overall = self.omnibus_test(name)
        pvalues = self.comparison_pvalues(name) if self.comparisons else {}
        for a, b in self.comparisons or []:
            if a not in set(self.sample_hap.haplotypes) or b not in set(self.sample_hap.haplotypes):
                raise ValueError('Comparison names must be present in the displayed haplotypes')

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
        jitter = np.random.RandomState(0).uniform(low=-0.25, high=0.25, size=len(ydata))
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
        span = max(float(np.ptp(ydata)), float(np.max(np.abs(ydata))) * 0.1, 1e-12)
        bracket_y = float(np.max(ydata))
        if self.comparisons is not None:
            for i, (group1, group2) in enumerate(self.comparisons):
                p_value = pvalues.get(tuple(sorted((group1, group2))), np.nan)
                sig_symbol = self.get_significance_symbol(p_value)
                if sig_symbol is None:
                    logger.warning(f'Skipping {group1} vs {group2}: comparison is not estimable (sample size or variance)')
                    continue
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
                    logger.info(f'{self.method}: {group1} vs {group2}, adjusted P={p_value:.17g}')
                    
                    # add significant line or bracket
                    max_y = bracket_y + 0.10 * span
                    if self.step_size is not None and i < len(self.step_size):
                        max_y += self.step_size[i]
                    bracket_y = max_y
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
        alignment = 'center' if rotation % 90 == 0 else ('right' if rotation > 0 else 'left')
        ax.set_xticklabels([f'{label}\n(n={len(data[label])})' for label in labels],
                           rotation=rotation, ha=alignment)
        if overall is not None:
            title = ax.set_title(overall['test'] + ': ' + self._p_label(overall['pvalue']), pad=12)
            title.set_gid('omnibus-test')
        ax.spines[["top", "right"]].set_visible(False)
        ax.yaxis.grid(
            True,
            linestyle='--',
            which='major',
            color='lightgrey',
            alpha=.5
        )
        save_figure(ax.figure, self.config['plot']['save_fig'], dpi=300, bbox_inches='tight')
