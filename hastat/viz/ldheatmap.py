import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.patches import RegularPolygon, Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable


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
    :param df: a square symmetric R^2 dataframe, with matching row/column labels
    :param plot_diag: whether to plot diagonal
    :param plot_value: whether to plot value
    :param plot_snp: whether to plot SNP location
    :param cmap: colormap, can be a string or a matplotlib colormap object
    :param ax: matplotlib axes object
    :return: matplotlib axes object
    """
    if not isinstance(df, pd.DataFrame) or df.shape[0] != df.shape[1]:
        raise ValueError('LD heatmap requires a square DataFrame')
    if df.empty:
        raise ValueError('No variants to plot in the LD heatmap')
    if not df.index.is_unique or not df.index.equals(df.columns):
        raise ValueError('LD heatmap row and column labels must match and be unique')
    data = df.to_numpy(dtype=float)
    if np.isinf(data).any() or np.any((data < 0) | (data > 1)):
        raise ValueError('LD values must lie in [0, 1] or be NaN')
    if not np.allclose(data, data.T, rtol=1e-6, atol=1e-8, equal_nan=True):
        raise ValueError('LD matrix must be symmetric, including missing values')
    if len(df) == 1 and not plot_diag:
        raise ValueError('A single variant has no pairwise cells; enable plot_diag')
    snp_loc = None
    if plot_snp and 'pos' in df.index.names:
        if 'chrom' in df.index.names and len(df.index.get_level_values('chrom').unique()) > 1:
            raise ValueError('A physical position track requires variants from one chromosome')
        snp_loc = df.index.get_level_values('pos').to_numpy(dtype=float)
        if not np.isfinite(snp_loc).all():
            raise ValueError('SNP positions must be finite')
    if ax is None:
        ax = plt.gca()

    # whether to plot diagonal
    start = 0
    stop = data.shape[0]
    n = data.shape[0]
    if not plot_diag:
        start = 1

    # set colormap
    cmap = mpl.colormaps.get_cmap(cmap).with_extremes(bad='#d1d5db')
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
    values = np.asarray(values, dtype=float)
    patch_collection.set_array(np.ma.masked_invalid(values))
    patch_collection.set_cmap(cmap)
    patch_collection.set_norm(norm)

    ax.add_collection(patch_collection)
    if np.isnan(values).any():
        ax.legend(handles=[Patch(facecolor='#d1d5db', label='Not estimable')],
                  loc='upper left', bbox_to_anchor=(0, -0.04), frameon=False, fontsize=8)
    ax.set_aspect('equal')
    ax.set_xlim(start * 0.5, stop - start * 0.5)
    ax.set_ylim(-0.1, (n - start) / 2 + 0.5 + 0.1)

    # Turn axis off
    ax.set_axis_off()

    # Add color bar
    cax = ax.inset_axes([1.03, 0.05, 0.035, 0.70])
    ax.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, shrink=.5, label=r"$R^2$")

    # Add text
    # Loop over data dimensions and create text annotations.
    # Change the text's color depending on the data.
    # Change the font's size depending on the patch.
    if plot_value:
        text_colors = ("black", "white")
        ax.figure.canvas.draw()
        cell_width_points = ax.get_window_extent().width * 72 / ax.figure.dpi / n
        font_size = max(1, min(10, cell_width_points / 4))
        for i, p in enumerate(patches):
            if np.isfinite(values[i]):
                ax.text(p.xy[0], p.xy[1], f'{values[i]:.2f}', ha='center', va='center',
                        fontsize=font_size, color=text_colors[int(values[i] > 0.5)])

    # Add SNP location
    # get SNP location from dataframe index name 'pos'
    if snp_loc is not None:
        span = np.ptp(snp_loc)
        mapped = (start * .5 + (snp_loc - snp_loc.min()) / span * (stop - start)
                  if span > 0 else np.full(n, (start * .5 + stop - start * .5) / 2))
        line_collection = LineCollection([[[x, 1], [i + 0.5, 0]] for i, x in enumerate(mapped)], linewidths=.5)
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
