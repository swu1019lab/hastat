"""Rendered-pixel regression checks for feature/node occlusion; run in PBS."""
import unittest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import test_gene_coordinates as coordinates
from hastat.viz.network import FancyNetwork


class PlotLayerValidation(unittest.TestCase):
    def test_gene_feature_centers_are_not_darkened_by_strand_line(self):
        for case in ['plus_flanks', 'minus_flanks']:
            fixture = coordinates.GeneCoordinateValidation()
            fixture.setUp()
            try:
                obj, _ = fixture.draw(case)
                pixels = np.asarray(obj.fig.canvas.buffer_rgba())
                for rect in [p for p in obj.ax.patches if isinstance(p, Rectangle)][:3]:
                    if rect.get_width() < .02:
                        continue
                    point = (rect.get_x()+rect.get_width()/2, rect.get_y()+rect.get_height()/2)
                    x, y = obj.ax.transData.transform(point)
                    observed = pixels[pixels.shape[0]-1-round(y), round(x), :3]
                    expected = np.rint(np.array(rect.get_facecolor()[:3])*255)
                    np.testing.assert_allclose(observed, expected, atol=1)
            finally:
                fixture.doCleanups()

    def test_translucent_nodes_hide_underlying_edges(self):
        obj = FancyNetwork.__new__(FancyNetwork)
        obj.group_colors = {'A': '#C5504B', 'B': '#90CAEE'}
        for composition in [{}, {'A': 1.0}, {'A': .65, 'B': .35}]:
            for background in ['white', '#e8eff4']:
                images = []
                for with_edge in [False, True]:
                    fig, ax = plt.subplots(figsize=(4, 4), dpi=100)
                    ax.set(xlim=(-1, 1), ylim=(-1, 1), facecolor=background)
                    if with_edge:
                        ax.plot([-1, 1], [0, 0], color='black', linewidth=5, zorder=1)
                    obj._draw_pie_node(ax, (0, 0), composition, 25, 500, .4)
                    fig.canvas.draw()
                    images.append(np.asarray(fig.canvas.buffer_rgba()).copy())
                    x, y = np.rint(ax.transData.transform((0, 0))).astype(int)
                    y = images[-1].shape[0]-1-y
                    plt.close(fig)
                np.testing.assert_array_equal(images[0][y-10:y+11, x-10:x+11],
                                              images[1][y-10:y+11, x-10:x+11])
                self.assertTrue(np.any(images[0] != images[1]))

