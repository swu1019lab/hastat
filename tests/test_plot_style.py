"""Palette, rendered legend geometry and visible strand-direction regressions."""
import tempfile
import unittest
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import networkx as nx
from hastat.viz.network import FancyNetwork
from hastat.viz.gene2 import FancyGene
from test_gene_coordinates import coordinate_inputs


class PlotStyleValidation(unittest.TestCase):
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
                scatter=ax.scatter([0],[0],s=[1])
                obj._create_legends(ax,scatter)
                fig.canvas.draw()
                renderer=fig.canvas.get_renderer()
                bounds=[c.get_window_extent(renderer) for c in obj.size_legend_markers]
                self.assertGreater(obj.legend_box.get_window_extent(renderer).x0,ax.bbox.x1)
                for count,box,label in zip(obj._size_legend_counts(),bounds,obj.size_legend_labels):
                    expected=np.sqrt(np.sqrt(count)*550)*fig.dpi/72
                    self.assertAlmostEqual(box.width,expected,places=6)
                    self.assertGreater(label.get_window_extent(renderer).x0,box.x1)
                for a,b in zip(bounds,bounds[1:]):
                    self.assertFalse(a.overlaps(b))
            finally:
                plt.close(fig)

    def test_strand_heads_are_visible_in_feature_gaps_for_both_strands(self):
        for case,sign in [('plus_flanks',1),('minus_flanks',-1)]:
            with tempfile.TemporaryDirectory() as tmp:
                cfg,_=coordinate_inputs(Path(tmp),case)
                obj=FancyGene(cfg); obj.plot(); obj.fig.canvas.draw()
                directions=[p for p in obj.ax.patches if p.get_gid()=='strand-direction']
                self.assertEqual(len(directions),2)
                for arrow in directions:
                    a,b=arrow._posA_posB
                    self.assertEqual(int(np.sign(b[0]-a[0])),sign)
                    self.assertGreater(abs(b[0]-a[0]),0)
                plt.close(obj.fig)

    def test_single_feature_direction_marker_is_horizontal_for_both_strands(self):
        for strand,sign in [('+',1),('-',-1)]:
            with tempfile.TemporaryDirectory() as tmp:
                cfg,_=coordinate_inputs(Path(tmp),'single_base')
                gff=Path(cfg['data']['gff_file'])
                gff.write_text(gff.read_text().replace('\t+\t',f'\t{strand}\t'))
                obj=FancyGene(cfg); obj.plot(); obj.fig.canvas.draw()
                marker=next(t for t in obj.ax.texts if t.get_gid()=='strand-direction')
                a,b=marker.arrow_patch._posA_posB
                self.assertEqual(int(np.sign(b[0]-a[0])),sign)
                self.assertAlmostEqual(a[1],b[1],places=6)
                self.assertGreater(abs(b[0]-a[0]),20)
                plt.close(obj.fig)

