import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.patches import RegularPolygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.transforms import Affine2D


def plot(
        df: pd.DataFrame,
        plot_diag: bool = True,
        plot_value: bool = False,
        plot_snp: bool = False,
        cmap: str = 'Reds',
        ax: plt.Axes = None):
    """
    LD heatmap plot

    Parameters
    ----------
    :param df: a dataframe with MxN R^2 values only
    :param plot_diag: whether to plot diagonal
    :param plot_value: whether to plot value
    :param plot_snp: whether to plot SNP location
    :param cmap: colormap, can be a string or a matplotlib colormap object
    :param ax: matplotlib axes object
    :return: matplotlib axes object
    """
    if ax is None:
        ax = plt.gca()

    # get R^2 values
    data = df.to_numpy()

    # whether to plot diagonal
    start = 0
    stop = data.shape[0]
    n = data.shape[0]
    if not plot_diag:
        start = 1

    # set colormap
    cmap = mpl.colormaps.get_cmap(cmap)
    norm = mpl.colors.Normalize(vmin=0, vmax=1)

    # get patch collection to plot
    patches = []
    values = []
    # for loop along the lower triangle of the matrix, then get the diagonal values
    for i in np.arange(start, stop):
        diag_values = np.diag(data, -i)
        values.extend(diag_values)
        for j in np.arange(0.5, len(diag_values) + 0.5):
            patches.append(RegularPolygon((j + i * 0.5, (n - i) / 2), numVertices=4, radius=0.5))

    patch_collection = PatchCollection(patches)
    patch_collection.set_array(values)
    patch_collection.set_cmap(cmap)
    patch_collection.set_norm(norm)

    ax.add_collection(patch_collection)
    ax.set_aspect('equal')
    ax.set_xlim(start * 0.5, stop - start * 0.5)
    ax.set_ylim(-0.1, (n - start) / 2 + 0.5 + 0.1)

    # Turn axis off
    ax.set_axis_off()

    # Add color bar
    cax = ax.inset_axes([0.8, 0.01, 0.03, 0.5])
    ax.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, shrink=.5, label=r"$R^2$")

    # Add text
    # Loop over data dimensions and create text annotations.
    # Change the text's color depending on the data.
    # Change the font's size depending on the patch.
    if plot_value:
        text_colors = ("black", "white")
        color_array = patch_collection.get_array()
        threshold = patch_collection.norm(color_array.max()) / 2
        for i, p in enumerate(patches):
            text = ax.text(p.xy[0], p.xy[1], "{:.2f}".format(values[i]),
                           ha="center", va="center",
                           color=text_colors[int(values[i] > threshold)])
            patch_bbox = p.get_window_extent()
            text_width = text.get_window_extent().transformed(ax.transData.inverted()).width
            font_size = text.get_fontsize()
            while font_size > 1 and text_width > patch_bbox.width / 2:
                font_size -= 1
                text.set_fontsize(font_size)
                text_width = text.get_window_extent().transformed(ax.transData.inverted()).width

    # Add SNP location
    # get SNP location from dataframe index name 'pos'
    if plot_snp and 'pos' in df.index.names:
        snp_loc = df.index.get_level_values('pos').to_numpy()
        sx = (stop - start) / (np.max(snp_loc) - np.min(snp_loc))
        scale_loc = Affine2D(). \
            translate(-np.min(snp_loc), ty=0). \
            scale(sx=sx, sy=1). \
            translate(tx=start * 0.5, ty=0). \
            transform(np.column_stack([snp_loc, [1] * n]))
        line_collection = LineCollection([[[a[0], 1], [i + 0.5, 0]] for i, a in enumerate(scale_loc)], linewidths=.5)
        line_collection.set_color('black')

        ax_divider = make_axes_locatable(ax)
        ax_snp = ax_divider.append_axes("top", size="10%", pad="0%", sharex=ax)
        ax_snp.add_collection(line_collection)
        ax_snp.set_xlim(start * 0.5, stop - start * 0.5)
        ax_snp.set_ylim(0, 1)
        ax_snp.set_xticks([])
        ax_snp.set_yticks([])
        ax_snp.set_yticklabels([])
        ax_snp.spines.bottom.set_visible(False)

        return ax_snp

    return ax
