# -*- coding: utf-8 -*-
# @Time    : 2024/10/26 12:58
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : gene.py


import pandas as pd
import numpy as np
from matplotlib import axes, patches, lines, text, ticker, collections
from collections import defaultdict
import matplotlib.pyplot as plt


class HapGene(object):
    def __init__(self, config: dict = None):
        if config is None:
            raise ValueError("config is None")
        self.config = config

        self.name = config['plot']['name']
        self.chrom = config['plot']['chrom']
        self.start = config['plot']['start']
        self.end = config['plot']['end']
        self.strand = config['plot']['strand']
        self.feature_data = dict()  # store feature data, e.g. exon, intron, utr, snp, gen, promoter
        self.feature_patches = defaultdict(list)
        self.feature_artists = defaultdict(list)
        self.feature_collections = defaultdict(list)
        self.legends = []
        self.center = 0.3  # center of gene

        # canvas visualize attributes
        self.x0 = config['plot']['start'] - (self.end - self.start) * 0.15
        self.x1 = config['plot']['end'] + (self.end - self.start) * 0.15
        self.y0 = 0
        self.y1 = 0.6

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
        _data = np.asarray(data)
        _data.sort(axis=0)
        self.add_feature_data('exon', _data)

    def add_introns(self, data: list or np.ndarray):
        """
        Add intron data to gene object

        :param data: A 2D list or numpy array contained intron data, format: [[start, end], [start, end], ...]
        :return:
        """
        _data = np.asarray(data)
        _data.sort(axis=0)
        self.add_feature_data('intron', _data)

    def add_utr(self, data: list or np.ndarray):
        """
        Add utr data to gene object

        :param data: A 2D list or numpy array contained utr data, format: [[start, end], [start, end], ...]
        :return:
        """
        _data = np.asarray(data)
        _data.sort(axis=0)
        self.add_feature_data('utr', _data)

    def add_snps(self, pos: list or np.ndarray, ref: list or np.ndarray, alt: list or np.ndarray,
                 gen: list or np.ndarray = None):
        """
        Add snp data to gene object

        :param pos: A 1D list or numpy array contained snp pos. format: [pos1, pos2, pos3, ...]
        :param ref: A 1D list or numpy array contained reference allele. format: 'A', 'T', 'C'
        :param alt: A 1D list or numpy array contained alternative allele. format: 'T', 'C', 'G', ...
        :param gen: A 2D list or numpy array contained genotype. format: [[0, 1, 2], [0, 1, 2], ...],
            (0: homozygous reference, 1: heterozygous, 2: homozygous alternative,
            each row is a sample, each column is a snp)
        :return:
        """
        index = np.argsort(pos)
        data = np.column_stack((np.asarray(pos)[index], np.asarray(ref)[index], np.asarray(alt)[index]))
        self.add_feature_data('snp', data)
        if gen is not None:
            self.add_feature_data('gen', np.asarray(gen)[:, index])

    def add_promoter(self, data: list or np.ndarray):
        """
        Add promoter data to gene object

        :param data: A 1D list or numpy array contained promoter data, format: [start, end]
        :return:
        """
        _data = np.asarray(data)
        _data.sort(axis=0)
        self.add_feature_data('promoter', _data)

    def draw_gene(self, fc='w', ec='C0', alpha=1, zorder=0, height=0.03):
        # define the center of gene with height and align parameters
        center = self.center
        patch = patches.Rectangle((self.start, center - height / 2), self.end - self.start, height, fc=fc, ec=ec,
                                  alpha=alpha, zorder=zorder)
        self.add_feature_patch('gene', patch)

    def draw_exons(self, fc='C0', ec='C0', alpha=1, zorder=4, height=0.03):
        if self.feature_data.get('exon') is None:
            raise ValueError("Exon data is not found")

        # define the center of exons with height and align parameters
        center = self.center
        for exon in self.feature_data['exon']:
            patch = patches.Rectangle((exon[0], center - height / 2), exon[1] - exon[0], height, fc=fc, ec=ec,
                                      alpha=alpha, zorder=zorder)
            self.add_feature_patch('exon', patch)

    def draw_introns(self, color='C5', alpha=1, zorder=2, shape='arrow', height=0.01):
        if self.feature_data.get('intron') is None:
            self.add_introns(self.get_introns_from_exons())

        # define the center of exons with height and align parameters
        center = self.center

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

    def draw_utr(self, fc='C7', ec='C7', alpha=1, zorder=3, height=0.03, align='center'):
        if self.feature_data.get('utr') is None:
            self.add_utr(self.get_utr_from_exons())

        # define the center of exons with height and align parameters
        center = self.center
        for utr in self.feature_data['utr']:
            patch = patches.Rectangle((utr[0], center - height / 2), utr[1] - utr[0], height, fc=fc, ec=ec,
                                      alpha=alpha, zorder=zorder)
            self.add_feature_patch('utr', patch)

    def draw_snps(self, offset=0.05, alpha=1, adjust=True, x0=None, x1=None, num=None):
        if adjust:
            if x0 is None:
                x0 = self.start
            if x1 is None:
                x1 = self.end
            if num is None:
                num = len(self.feature_data['snp']) + 1
            intervals = np.linspace(x0, x1, num)
            # get the center of each interval
            new_pos = (intervals[:-1] + intervals[1:]) / 2
        else:
            new_pos = self.feature_data['snp'][:, 0]

        center = self.center
        for i, (pos, ref, alt) in enumerate(self.feature_data['snp']):
            ann1 = text.Annotation(ref, (float(new_pos[i]), self.y1 + offset), xycoords="data",
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

    def draw_genotype(self, x0=None, x1=None, y0=-1, y1=0, alpha=1, zorder=5,
                      row_grids: list = None, row_labels: dict = None):
        """
        Draw genotype heatmap

        :param x0: start position of genotype heatmap in x-axis
        :param x1: end position of genotype heatmap in x-axis
        :param y0: start position of genotype heatmap in y-axis
        :param y1: end position of genotype heatmap in y-axis
        :param alpha: transparency of genotype heatmap
        :param zorder: zorder of genotype heatmap
        :param row_grids: a list contained row grids
        :param row_labels: a dict contained row index and labels, e.g. {0: 'A', 1: 'B'}
        :return:
        """
        rows, cols = self.feature_data['gen'].shape
        snp_pos = self.feature_data['snp'][:, 0].astype(float)
        center = self.center
        x0 = min(self.x0, np.min(snp_pos)) if x0 is None else x0
        x1 = max(self.x1, np.max(snp_pos)) if x1 is None else x1
        x = np.linspace(x0, x1, cols + 1).tolist()
        y = np.linspace(y0, y1, rows + 1).tolist()[::-1]
        w = x[1] - x[0]
        h = y[1] - y[0]
        colors = {0: 'steelblue', 1: 'tan', 2: 'tomato'}
        rects = []
        # draw heatmap
        for i in range(rows):
            for j in range(cols):
                color = colors.get(self.feature_data['gen'][i, j], 'w')
                rect = patches.Rectangle((x[j], y[i]), w, h, fc=color, ec=color, alpha=alpha, zorder=zorder)
                rects.append(rect)
        self.add_feature_collections('gen', collections.PatchCollection(rects, match_original=True))
        # draw snp line
        for i in range(cols):
            line = lines.Line2D([snp_pos[i], x[i] + w / 2], [center, y[0]], color='C7',
                                zorder=zorder, alpha=alpha, linestyle='--')
            self.add_feature_artists('link', line)
        # draw grid and labels
        if row_grids is not None and row_labels is not None and len(row_grids) == len(row_labels):
            for row in row_grids:
                line = lines.Line2D([x[0], x[-1]], [y[row], y[row]],
                                    color='w', zorder=zorder + 1, alpha=alpha)
                self.add_feature_artists('row_grids', line)
            for row, label in row_labels.items():
                ann = text.Annotation(label, (x[0], y[int(row)]), xycoords="data",
                                      xytext=(-10, 0), textcoords="offset points",
                                      va="center", ha="right", color='w',
                                      bbox=dict(boxstyle="round", fc="C1", ec="w", alpha=alpha), zorder=zorder + 1)
                # txt = text.Text(x[0], y[int(row)], label, ha='right', va='center', color='k', zorder=zorder + 1)
                self.add_feature_artists('row_labels', ann)
        # recalculate the size of canvas
        self.x0 = min(self.x0, np.min(snp_pos))
        self.x1 = max(self.x1, np.max(snp_pos))
        self.y0 = min(self.y0, y0)
        self.y1 = max(self.y1, y1)

    def draw_promoter(self, fc='C7', ec='C7', alpha=1, zorder=1, height=0.01):
        if self.feature_data.get('promoter') is None or len(self.feature_data['promoter']) == 0:
            self.add_promoter(self.get_promoter())

        center = self.center
        promoter = self.feature_data['promoter']
        if np.asarray(promoter).ndim == 2:
            promoter = promoter[0]
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

    def add_feature_collections(self, feature: str, collection):
        self.feature_collections[feature] = collection

    def update_feature_patch(self, feature: str, index: int, patch: patches.Patch):
        self.feature_patches[feature][index] = patch

    def add_data(self):
        df1 = pd.read_csv(self.config['data']['ann_file'], header=0)
        # 0-chrom, 1-start, 2-end, 3-name, 4-strand, 5-type
        self.add_exons(df1.query('type == "exon" & name == @self.name').iloc[:, [1, 2]].to_numpy(dtype=np.int64))
        self.add_introns(df1.query('type == "intron" & name == @self.name').iloc[:, [1, 2]].to_numpy(dtype=np.int64))
        self.add_utr(df1.query('type == "utr" & name == @self.name').iloc[:, [1, 2]].to_numpy(dtype=np.int64))
        self.add_promoter(df1.query('type == "promoter" & name == @self.name').iloc[:, [1, 2]].to_numpy(dtype=np.int64))

        df2 = pd.read_csv(self.config['data']['hap_file'], header=[0, 1, 2, 3], index_col=0)
        df2 = df2[df2.loc[:, ('haplotypes', slice(None), slice(None), slice(None))].iloc[:, 0].isin(
            self.config['plot']['show_hap'])]
        df2 = df2.sort_values(by=df2.columns[-1], axis=0)
        hap, count = np.unique(df2.iloc[:, -1], return_counts=True)
        self.add_snps(
            pos=df2.columns.get_level_values('pos')[:-1].to_numpy(dtype=np.int64),
            ref=df2.columns.get_level_values('ref')[:-1].tolist(),
            alt=df2.columns.get_level_values('alt')[:-1].tolist(),
            gen=df2.iloc[:, :-1].values
        )
        self.add_feature_data('row_labels', dict(zip(np.cumsum(count) - count / 2, hap)))
        self.add_feature_data('row_grids', np.cumsum(count))

    def plot(self, ax: axes.Axes = None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(self.config['plot']['width'],
                                            self.config['plot']['height']))

        self.add_data()

        self.draw_exons()
        self.draw_introns()
        self.draw_utr()
        self.draw_snps()
        self.draw_genotype(
            row_grids=self.feature_data['row_grids'],
            row_labels=self.feature_data['row_labels']
        )
        self.draw_promoter()

        ax.autoscale(enable=True, axis='both')
        for feature in self.feature_patches:
            for patch in self.feature_patches[feature]:
                ax.add_patch(patch)
        for feature in self.feature_artists:
            for artist in self.feature_artists[feature]:
                ax.add_artist(artist)
        for feature in self.feature_collections:
            ax.add_collection(self.feature_collections[feature])

        ax.set_ylim(self.y0 - 0.15, self.y1 + 0.15)
        ax.set_xlim(self.x0 - 2000, self.x1 + 2000)

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
        plt.savefig(self.config['plot']['save_fig'], dpi=300)
