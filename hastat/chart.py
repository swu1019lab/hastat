# -*- coding: utf-8 -*-
# @Time    : 2024/5/9 11:47
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : chart.py

import numpy as np
import pandas as pd
from matplotlib import axes, patches, lines, text, ticker
from collections import defaultdict, Counter
from cycler import cycle


class HapBar(object):
    def __init__(self):
        self.hap_data = None
        self.group_data = None
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
            self.hap_data = df
        else:
            self.hap_data = df[df['haplotypes'].isin(labels)]

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
        self.group_data = df

    def merge_data(self):
        if self.hap_data is None or self.group_data is None:
            raise ValueError("hap_data or group_data is None")
        self.data = self.hap_data.merge(self.group_data, how='inner', on='sample')

    def get_plot_data(self, name: str, calc_percentage: bool = False):
        data = self.data.loc[:, ['sample', 'haplotypes', name]].groupby(['haplotypes', name]).size().unstack(
            fill_value=0)
        if calc_percentage:
            return data.div(data.sum())
        return data

    def plot(self, ax: axes.Axes, name: str = None):
        if ax is None:
            raise ValueError("ax is None")
        self.merge_data()
        if name is None:
            name = self.group_data.columns[1]
        colors = cycle(['#498DCB', '#F9BEBF', '#747474', '#EE3424'])
        data = self.get_plot_data(name, calc_percentage=True)
        bottom = np.zeros(data.shape[1])
        x = data.columns.get_level_values(name).to_numpy()
        for i, (hap, row) in enumerate(data.iterrows()):
            ax.bar(x, row.to_numpy(), width=0.5, bottom=bottom, label=hap, color=next(colors))
            bottom += row.to_numpy()
        ax.set_ylabel('Haplotype frequency')
        ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        ax.yaxis.grid(
            True,
            linestyle='--',
            which='major',
            color='lightgrey',
            alpha=.5
        )


class HapBox(object):
    def __init__(self):
        self.hap_data = None
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
            self.hap_data = df
        else:
            self.hap_data = df[df['haplotypes'].isin(labels)]

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
        if self.hap_data is None or self.phe_data is None:
            raise ValueError("hap_data or phe_data is None")
        self.data = self.hap_data.merge(self.phe_data, how='inner', on='sample')

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


class Gene(object):
    def __init__(self, name: str, chrom: str, start: float, end: float, strand: str):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.feature_data = dict()
        self.feature_patches = defaultdict(list)
        self.feature_artists = defaultdict(list)
        self.legends = []

        # visualize attributes
        self.y0 = 0
        self.y1 = 0.4

    def __str__(self):
        return f"Gene(name={self.name})"

    def __repr__(self):
        return f"Gene(name={self.name})"

    def __len__(self):
        return self.end - self.start + 1

    def __getitem__(self, key):
        return self.feature_data[key]

    def add_exons(self, data: list or np.ndarray):
        """
        Add exon data to gene object

        :param data: A 2D list or numpy array contained exon data, format: [[start, end], [start, end], ...]
        :return:
        """
        data = np.asarray(data)
        data.sort(axis=0)
        self.add_feature_data('exon', data)

    def add_introns(self, data: list or np.ndarray):
        """
        Add intron data to gene object

        :param data: A 2D list or numpy array contained intron data, format: [[start, end], [start, end], ...]
        :return:
        """
        data = np.asarray(data)
        data.sort(axis=0)
        self.add_feature_data('intron', data)

    def add_utr(self, data: list or np.ndarray):
        """
        Add utr data to gene object

        :param data: A 2D list or numpy array contained utr data, format: [[start, end], [start, end], ...]
        :return:
        """
        data = np.asarray(data)
        data.sort(axis=0)
        self.add_feature_data('utr', data)

    def add_snps(self, pos: list or np.ndarray, ref: list or np.ndarray, alt: list or np.ndarray):
        """
        Add snp data to gene object

        :param pos: A 1D list or numpy array contained snp pos. format: [pos1, pos2, pos3, ...]
        :param ref: A 1D list or numpy array contained reference allele. format: 'A', 'T', 'C'
        :param alt: A 1D list or numpy array contained alternative allele. format: 'T', 'C', 'G', ...
        :return:
        """
        index = np.argsort(pos)
        data = np.column_stack((np.asarray(pos)[index], np.asarray(ref)[index], np.asarray(alt)[index]))
        self.add_feature_data('snp', data)

    def add_promoter(self, data: list or np.ndarray):
        """
        Add promoter data to gene object

        :param data: A 1D list or numpy array contained promoter data, format: [start, end]
        :return:
        """
        data = np.asarray(data)
        data.sort(axis=0)
        self.add_feature_data('promoter', data)

    def draw_exons(self, fc='C0', ec='C0', alpha=1, zorder=4, height=0.03):
        center = (self.y0 + self.y1) / 2
        for exon in self.feature_data['exon']:
            patch = patches.Rectangle((exon[0], center - height / 2), exon[1] - exon[0], height, fc=fc, ec=ec,
                                      alpha=alpha, zorder=zorder)
            self.add_feature_patch('exon', patch)

    def draw_introns(self, color='C5', alpha=1, zorder=2, shape='polyline', height=0.01):
        if self.feature_data.get('intron') is None:
            self.add_introns(self.get_introns_from_exons())

        center = (self.y0 + self.y1) / 2
        for intron in self.feature_data['intron']:
            if shape == 'polyline':
                line0 = lines.Line2D([intron[0], (intron[0] + intron[1]) / 2],
                                     [center, center + height / 2],
                                     color=color, zorder=zorder, alpha=alpha)

                line1 = lines.Line2D([(intron[0] + intron[1]) / 2, intron[1]],
                                     [center + height / 2, center],
                                     color=color, zorder=zorder, alpha=alpha)

                self.add_feature_artists('intron', line0)
                self.add_feature_artists('intron', line1)

            elif shape == 'arrow':
                if self.strand == '+':
                    self.add_feature_patch('intron', patches.FancyArrowPatch((intron[0], center),
                                                                             (intron[1], center),
                                                                             color=color, zorder=zorder, alpha=alpha,
                                                                             arrowstyle='->', mutation_scale=10,
                                                                             shrinkA=0, shrinkB=0))
                else:
                    self.add_feature_patch('intron', patches.FancyArrowPatch((intron[1], center),
                                                                             (intron[0], center),
                                                                             color=color, zorder=zorder, alpha=alpha,
                                                                             arrowstyle='<-', mutation_scale=10,
                                                                             shrinkA=0, shrinkB=0))
            else:
                line = lines.Line2D([intron[0], intron[1]],
                                    [center, center],
                                    color=color, zorder=zorder, alpha=alpha, linestyle='-',
                                    solid_capstyle='butt')
                self.add_feature_artists('intron', line)

    def draw_utr(self, fc='C7', ec='C7', alpha=1, zorder=3, height=0.03):
        if self.feature_data.get('utr') is None:
            self.add_utr(self.get_utr_from_exons())

        center = (self.y0 + self.y1) / 2
        for utr in self.feature_data['utr']:
            patch = patches.Rectangle((utr[0], center - height / 2), utr[1] - utr[0], height, fc=fc, ec=ec,
                                      alpha=alpha, zorder=zorder)
            self.add_feature_patch('utr', patch)

    def draw_snps(self, height=0.15, alpha=1, adjust=True, start=None, end=None, num=None):
        if adjust:
            if start is None:
                start = self.start
            if end is None:
                end = self.end
            if num is None:
                num = len(self.feature_data['snp']) + 1
            intervals = np.linspace(start, end, num)
            # get the center of each interval
            new_pos = (intervals[:-1] + intervals[1:]) / 2
        else:
            new_pos = self.feature_data['snp'][:, 0]

        center = (self.y0 + self.y1) / 2
        for i, (pos, ref, alt) in enumerate(self.feature_data['snp']):
            ann1 = text.Annotation(ref, (float(new_pos[i]), center + height), xycoords="data",
                                   va="center", ha="center", color='w',
                                   bbox=dict(boxstyle="round", fc='C0', ec="C0", alpha=alpha))
            offset_from = text.OffsetFrom(ann1, (0.5, 0))
            ann2 = text.Annotation(alt, (float(pos), center), xycoords="data",
                                   xytext=(0, -10), textcoords=offset_from,
                                   va="top", ha="center", color='w',
                                   bbox=dict(boxstyle="round", fc="C3", ec="C3", alpha=alpha),
                                   arrowprops=dict(arrowstyle="-", color="C7", alpha=alpha))
            self.add_feature_artists('snp', ann1)
            self.add_feature_artists('snp', ann2)

    def draw_promoter(self, fc='C7', ec='C7', alpha=1, zorder=1, height=0.01):
        if self.feature_data.get('promoter') is None:
            self.add_promoter(self.get_promoter())

        center = (self.y0 + self.y1) / 2
        promoter = self.feature_data['promoter']
        patch = patches.Rectangle((promoter[0], center - height / 2), promoter[1] - promoter[0], height,
                                  fc=fc, ec=ec, alpha=alpha, zorder=zorder)
        self.add_feature_patch('promoter', patch)

    def get_exons(self):
        return self.feature_data['exon']

    def get_promoter(self, length=2000):
        if self.strand == '+':
            return [self.start - length, self.start]
        else:
            return [self.end, self.end + length]

    def get_utr_from_exons(self):
        """
        Get utr location from exons

        :return: utr array
        """
        utr = [[self.start, self.feature_data['exon'][0][0]], [self.feature_data['exon'][-1][1], self.end]]
        return np.asarray(utr)

    def get_introns_from_exons(self):
        """
        Get introns location of gene

        :return: intron array
        """
        starts = self.feature_data['exon'][:-1, 1]
        ends = self.feature_data['exon'][1:, 0]
        return np.column_stack((starts, ends))

    def reverse_strand(self):
        """
        Reverse the direction of strand (from 3'->5' to 5'->3') to get new position

        :return: None
        """
        if self.strand == '-':
            self.feature_data['exon'] = self.end - self.feature_data['exon'][::-1][:, ::-1]
            self.feature_data['exon'] += self.start

    def set_zero_start(self):
        """
        Set start position of gene to zero

        :return: None
        """
        self.set_offset(-self.start, 0)

    def set_offset(self, offset_x=None, offset_y=None):
        """
        Set offset of gene

        :param offset_x: offset value in x-axis
        :param offset_y: offset value in y-axis
        :return: None
        """
        self.start += offset_x
        self.end += offset_x
        self.feature_data['exon'] += offset_x
        self.feature_data['utr'] += offset_x
        self.y0 += offset_y
        self.y1 += offset_y

    def add_feature_data(self, feature: str, data: np.ndarray or list):
        self.feature_data[feature] = data

    def add_feature_patch(self, feature: str, patch: patches.Patch):
        self.feature_patches[feature].append(patch)

    def add_feature_artists(self, feature: str, artist):
        self.feature_artists[feature].append(artist)

    def update_feature_patch(self, feature: str, index: int, patch: patches.Patch):
        self.feature_patches[feature][index] = patch

    def plot(self, ax: axes.Axes):
        if ax is None:
            raise ValueError("ax is None")
        ax.autoscale(enable=True, axis='both')
        for feature in self.feature_patches:
            for patch in self.feature_patches[feature]:
                ax.add_patch(patch)
        for feature in self.feature_artists:
            for artist in self.feature_artists[feature]:
                ax.add_artist(artist)
        ax.set_ylim(self.y0, self.y1)

        # set x-axis ticks with kb unit
        ax.set_title('Position (Mb)')
        ax.xaxis.set_major_formatter(
            ticker.FuncFormatter(lambda x, pos: '{:.3f}'.format(x / 1000000))
        )
        ax.spines[["left", "bottom", "right"]].set_visible(False)
        ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, left=False, labelleft=False)
        # set grid
        ax.xaxis.grid(
            True,
            linestyle='--',
            which='major',
            color='lightgrey',
            alpha=.5
        )
        # Hide these grid behind plot objects
        ax.set_axisbelow(True)
