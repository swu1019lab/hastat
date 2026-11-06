"""Palette, rendered legend geometry and visible strand-direction regressions."""
import tempfile
import unittest
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd
import networkx as nx
from hastat.viz.network import FancyNetwork
from hastat.viz.gene import FancyGene
from test_gene_coordinates import coordinate_inputs


class PlotStyleValidation(unittest.TestCase):
    def draw_network_edge(self, endpoint, symbol='|', weight=6, style='line', show=True):
        obj = FancyNetwork.__new__(FancyNetwork)
        obj.network = nx.Graph()
        obj.network.add_edge('A', 'B', weight=weight)
        obj.config = {'plot': {'show_edge_label': show, 'edge_label_symbol': symbol,
                               'edge_label_size': 8, 'edge_style': style}}
        obj.node_sizes = {'A': 1, 'B': 1}
        obj.hap_composition = {}
        obj.group_colors = {}
        fig = plt.figure(figsize=(5, 3), dpi=100)
        self.addCleanup(plt.close, fig)
        ax = fig.add_axes([0, 0, 1, 1], xlim=(0, 200), ylim=(0, 200))
        positions = {'A': np.array([30., 30.]), 'B': np.array(endpoint, dtype=float)}
        obj._draw_graph(ax, positions, {'A': 8., 'B': 8.}, 256)
        fig.canvas.draw()
        marks = [c for c in ax.collections if isinstance(c, LineCollection)]
        return fig, ax, marks

    def test_edge_ticks_count_position_and_perpendicular_direction(self):
        for endpoint in ([170, 30], [30, 170], [170, 150]):
            for style in ('line', 'curve'):
                with self.subTest(endpoint=endpoint, style=style):
                    fig, ax, marks = self.draw_network_edge(endpoint, style=style)
                    self.assertEqual(len(marks), 1)
                    self.assertFalse(ax.texts)
                    segments = np.array(marks[0].get_segments())
                    self.assertEqual(len(segments), 6)
                    display = ax.transData.transform(segments.reshape(-1, 2)).reshape(-1, 2, 2)
                    path = ax.transData.transform(ax.lines[0].get_xydata())
                    vectors = np.diff(path, axis=0)
                    centers = display.mean(axis=1)
                    lengths = np.linalg.norm(vectors, axis=1)
                    arcs = []
                    for center, segment in zip(centers, display):
                        fractions = np.clip(np.sum((center-path[:-1])*vectors, axis=1) / lengths**2, 0, 1)
                        projections = path[:-1] + fractions[:, None]*vectors
                        nearest = np.argmin(np.linalg.norm(projections-center, axis=1))
                        np.testing.assert_allclose(projections[nearest], center, atol=1e-8)
                        self.assertAlmostEqual(np.dot(segment[1]-segment[0], vectors[nearest]), 0, places=7)
                        arcs.append(lengths[:nearest].sum()+fractions[nearest]*lengths[nearest])
                    np.testing.assert_allclose(np.diff(arcs), np.diff(arcs)[0], atol=1e-8)
                    self.assertGreater(arcs[0], 0)
                    self.assertLess(arcs[-1], lengths.sum())
                    self.assertAlmostEqual(arcs[0]+arcs[-1], lengths.sum(), places=7)

    def test_edge_ticks_keep_weights_and_label_visibility(self):
        for weight, expected in [(0, 0), (1, 1), (18, 18)]:
            with self.subTest(weight=weight):
                _, ax, marks = self.draw_network_edge([170, 150], weight=weight)
                self.assertEqual(sum(len(c.get_segments()) for c in marks), expected)
                self.assertFalse(ax.texts)
        _, ax, marks = self.draw_network_edge([170, 150], show=False)
        self.assertFalse(marks)
        self.assertFalse(ax.texts)
        for symbol, weight, label in [(None, 6, '6'), ('*', 3, '***'), ('|', 2.5, '2.5')]:
            with self.subTest(symbol=symbol, weight=weight):
                _, ax, marks = self.draw_network_edge([170, 150], symbol=symbol, weight=weight)
                self.assertFalse(marks)
                self.assertEqual([t.get_text() for t in ax.texts], [label])

    def test_default_network_palette_matches_pi_group_order(self):
        obj=FancyNetwork.__new__(FancyNetwork)
        groups=['Panel1','Panel2','Panel3']
        colors=obj._generate_colors(groups)
        expected=LinearSegmentedColormap.from_list('pi_cmap',
            ['#C5504B','#114F8B','#FCE988','#90CAEE'],N=100)(np.linspace(0,1,3))[:,:3]
        np.testing.assert_allclose([colors[g] for g in groups],expected)
        custom=obj._generate_colors(groups,['navy'])
        np.testing.assert_allclose(custom['Panel1'],[0,0,128/255])
        np.testing.assert_allclose(custom['Panel2'],expected[1])

    def make_network(self, counts):
        obj=FancyNetwork.__new__(FancyNetwork)
        obj.network=nx.Graph()
        obj.network.add_nodes_from(range(len(counts)))
        obj.node_sizes=dict(enumerate(counts))
        obj.sample_hap_data=pd.DataFrame({'samples':['demo']})
        obj.group_colors={'Panel1':'#C5504B','Panel2':'#90CAEE'}
        obj.config={'plot':{'node_scale_factor':550,'node_alpha':1}}
        return obj

    def test_size_legend_uses_actual_counts_and_omits_unknown_counts(self):
        obj=self.make_network([36,30,24,20,18,14,10,8,0])
        self.assertEqual(obj._size_legend_counts(),[8,20,36])
        obj.node_sizes=dict.fromkeys(obj.network,7)
        self.assertEqual(obj._size_legend_counts(),[7])
        obj.node_sizes=dict.fromkeys(obj.network,0)
        self.assertEqual(obj._size_legend_counts(),[])
        obj.sample_hap_data=None
        self.assertEqual(obj._size_legend_counts(),[])

    def test_size_legend_circles_match_nodes_and_never_touch_labels(self):
        for counts in [[36,30,24,20,18,14,10,8],[5],[1,2,3]]:
            obj=self.make_network(counts)
            fig,ax=plt.subplots(figsize=(5,4),dpi=100)
            try:
                obj._legend(ax)
                fig.canvas.draw()
                renderer=fig.canvas.get_renderer()
                bounds=[c.get_window_extent(renderer) for c in obj.size_legend_markers]
                for count,box,label in zip(obj._size_legend_counts(),bounds,obj.size_legend_labels):
                    expected=np.sqrt(count*550)*fig.dpi/72
                    self.assertAlmostEqual(box.width,expected,places=6)
                    self.assertGreater(label.get_window_extent(renderer).x0,box.x1)
                for a,b in zip(bounds,bounds[1:]):
                    self.assertAlmostEqual(a.y0,b.y0,places=6)
                    self.assertAlmostEqual((a.x0+a.x1)/2,(b.x0+b.x1)/2,places=6)
            finally:
                plt.close(fig)

    def test_single_terminal_marker_for_both_strands(self):
        for case,sign in [('plus_flanks',1),('minus_flanks',-1)]:
            with tempfile.TemporaryDirectory() as tmp:
                cfg,_=coordinate_inputs(Path(tmp),case)
                obj=FancyGene(cfg); obj.plot(); obj.fig.canvas.draw()
                markers=[t for t in obj.ax.texts if t.get_gid()=='strand-direction']
                self.assertEqual(len(markers),1)
                self.assertFalse([p for p in obj.ax.patches if p.get_gid()=='strand-direction'])
                marker=markers[0]
                self.assertEqual(marker.get_text(), '>' if sign == 1 else '<')
                self.assertIsNone(marker.arrow_patch)
                anchor=obj.ax.transData.transform(marker.xy)
                blocks=[p.get_window_extent() for p in obj.ax.patches
                        if isinstance(p,matplotlib.patches.Rectangle)][:len(obj.gene_data['CDS'])]
                edge=max(p.xmax for p in blocks) if sign==1 else min(p.xmin for p in blocks)
                self.assertAlmostEqual(anchor[0],edge,places=6)
                bounds=marker.get_window_extent()
                gap=bounds.xmin-edge if sign==1 else edge-bounds.xmax
                self.assertGreater(gap,0)
                self.assertEqual(marker.get_position(),(3*sign,0))
                plt.close(obj.fig)

    def test_single_feature_direction_marker_is_horizontal_for_both_strands(self):
        for strand,sign in [('+',1),('-',-1)]:
            with tempfile.TemporaryDirectory() as tmp:
                cfg,_=coordinate_inputs(Path(tmp),'single_base')
                gff=Path(cfg['data']['gff_file'])
                gff.write_text(gff.read_text().replace('\t+\t',f'\t{strand}\t'))
                obj=FancyGene(cfg); obj.plot(); obj.fig.canvas.draw()
                markers=[t for t in obj.ax.texts if t.get_gid()=='strand-direction']
                self.assertEqual(len(markers),1)
                self.assertEqual(markers[0].get_text(), '>' if sign==1 else '<')
                self.assertIsNone(markers[0].arrow_patch)
                self.assertEqual(markers[0].get_position(),(3*sign,0))
                plt.close(obj.fig)
