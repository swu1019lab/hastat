"""Independent inclusive-coordinate and rendered-pixel checks; execute in PBS."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
import numpy as np
import pandas as pd

from hastat.viz.gene2 import FancyGene


def coordinate_inputs(folder, case='plus_flanks'):
    folder.mkdir(parents=True, exist_ok=True)
    strand = '-' if case == 'minus_flanks' else '+'
    single = case == 'single_base'
    lo, hi = (100, 100) if single else (100, 200)
    upstream, downstream = (0, 0) if single else (50, 20)
    start, end = ((lo-downstream, hi+upstream) if strand == '-' else (lo-upstream, hi+downstream))
    positions = [100] if single else [start, 100, 150, 200, end]
    (folder/'gene.gff3').write_text(
        f'1\ttest\tgene\t{lo}\t{hi}\t.\t{strand}\t.\tID=G1\n' +
        ''.join(f'1\ttest\tCDS\t{a}\t{b}\t.\t{strand}\t0\tID=C{i};Parent=G1\n'
                for i,(a,b) in enumerate([(100,100)] if single else [(100,120),(150,150),(180,200)])))
    table = pd.DataFrame(dict(chrom=['1']*len(positions),pos=positions,ref='A',alt='G',
                             Hap1=[0]*len(positions),Hap2=[2]*len(positions)))
    table.to_csv(folder/'gene.table.csv',index=False)
    (folder/'gene.group.csv').write_text('samples,haplotypes\ns1,Hap1\ns2,Hap1\ns3,Hap2\n')
    pd.DataFrame(dict(chrom='1',pos=positions,ann=[f'v{p}' for p in positions])).to_csv(folder/'ann.csv',index=False)
    for stat, values in [('pi',[.001] if single else [.001,.002,np.nan,.004,.005]),
                         ('fst',[-.1] if single else [-.1,.1,np.nan,.2,.3])]:
        pd.DataFrame(dict(chrom='1',start=positions,end=positions,**{stat:values})).to_csv(folder/f'{stat}.csv',index=False)
    custom = case == 'custom_frame'
    (folder/'gene.toml').write_text('[default]\nbbox=[0.1,0.1,0.8,0.8]\n' +
        ('xy=[0.15,0.78]\nwidth=0.65\n' if custom else 'xy=[0.0,0.78]\nwidth=1.0\n') +
        'height=0.035\nann_y=0.92\nhap_y=0.52\nhap_height=0.12\nhap_freq_show=false\n'
        'hap_file="gene"\nann_file="ann.csv"\npi_file="pi.csv"\nfst_file="fst.csv"\n'
        'pi_y=0.26\npi_height=0.14\nfst_y=0.0\nfst_height=0.14\n')
    config = dict(data=dict(gff_file=str(folder/'gene.gff3'),genes=['G1'],toml_file=str(folder/'gene.toml')),
                  plot=dict(width=10,height=7,upstream=upstream,downstream=downstream))
    return config, table


class GeneCoordinateValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.addCleanup(plt.close, 'all')
        self.root = Path(self.tmp.name)

    def draw(self, case='plus_flanks'):
        cfg, table = coordinate_inputs(self.root/case,case)
        obj = FancyGene(cfg)
        obj.plot()
        obj.fig.canvas.draw()
        return obj, table

    def test_inclusive_feature_widths_and_single_base_feature(self):
        obj, _ = self.draw()
        rectangles = [p for p in obj.ax.patches if isinstance(p,Rectangle)][:3]
        np.testing.assert_allclose([p.get_width() for p in rectangles], .8*np.array([21,1,21])/171)
        np.testing.assert_allclose([p.get_x() for p in rectangles], .1+.8*np.array([50,100,130])/171)

    def test_plus_and_minus_strand_flanks_and_arrow_direction(self):
        for case, interval, direction in [('plus_flanks',(50,220),1),('minus_flanks',(80,250),-1)]:
            obj,_ = self.draw(case)
            gene = obj.gene_data['G1']
            self.assertEqual((gene['start'],gene['end']),interval)
            arrow = next(p for p in obj.ax.patches if isinstance(p,FancyArrowPatch))
            self.assertEqual(np.sign(arrow._posA_posB[1][0]-arrow._posA_posB[0][0]),direction)
            label = next(t for t in obj.ax.texts if t.get_text() == 'G1')
            self.assertAlmostEqual(label.xy[0],.1+.8*(150-interval[0]+.5)/171)
            self.assertLess(*obj.ax.child_axes[0].get_xlim())

    def test_rendered_pixel_alignment_across_all_physical_tracks(self):
        for case in ['plus_flanks','minus_flanks','custom_frame','single_base']:
            obj,table = self.draw(case)
            start,end = (obj.gene_data['G1'][k] for k in ['start','end'])
            x0,width = (.22,.52) if case == 'custom_frame' else (.1,.8)
            for i,pos in enumerate(table.pos):
                expected_x = x0+width*(pos-start+.5)/(end-start+1)
                connector = obj.ax.lines[2*i+1]
                np.testing.assert_allclose(connector.get_xdata(),[expected_x,expected_x],atol=1e-12)
                annotation = next(t for t in obj.ax.texts if t.get_text() == f'v{pos}')
                self.assertAlmostEqual(annotation.xy[0],expected_x)
                expected_px = obj.ax.transData.transform((expected_x,0))[0]
                for track in obj.ax.child_axes:
                    self.assertAlmostEqual(track.transData.transform((pos,0))[0],expected_px,places=8)

    def test_haplotype_columns_remain_equal_spaced_and_counts_unchanged(self):
        obj, table = self.draw('custom_frame')
        np.testing.assert_allclose([line.get_xdata()[0] for line in obj.ax.lines[::2]],
                                  .1+.8*(np.arange(len(table))+.5)/len(table))
        self.assertEqual([t.get_text() for t in obj.ax.texts if t.get_text().startswith('Hap')],
                         ['Hap1 (n=2)','Hap2 (n=1)'])

    def test_connector_reaches_gene_lower_edge(self):
        obj,_ = self.draw()
        bottom = .1+.8*(.78-.035/2)
        for diagonal,vertical in zip(obj.ax.lines[::2],obj.ax.lines[1::2]):
            self.assertAlmostEqual(vertical.get_ydata()[-1],bottom)
            self.assertAlmostEqual(diagonal.get_ydata()[-1],(bottom+(.1+.8*(.52+.12)))/2)

    def test_single_base_gene_full_plot_no_zero_division(self):
        obj,_ = self.draw('single_base')
        rect = next(p for p in obj.ax.patches if isinstance(p,Rectangle))
        self.assertAlmostEqual(rect.get_width(),.8)
        for track in obj.ax.child_axes:
            self.assertEqual(track.get_xlim(),(99.5,100.5))
        self.assertAlmostEqual(obj.ax.lines[1].get_xdata()[0],.5)

    def test_flank_annotations_retained_but_outside_region_excluded(self):
        cfg,_ = coordinate_inputs(self.root)
        with (self.root/'ann.csv').open('a') as stream:
            stream.write('1,49,outside_left\n1,221,outside_right\n2,100,other_chrom\n')
        obj = FancyGene(cfg)
        obj.plot()
        labels = [t.get_text() for t in obj.ax.texts]
        self.assertIn('v50',labels)
        self.assertIn('v220',labels)
        self.assertFalse(any(label in labels for label in ['outside_left','outside_right','other_chrom']))

    def test_window_values_nan_gaps_and_negative_fst_unchanged(self):
        obj,table = self.draw('custom_frame')
        for track, expected in zip(obj.ax.child_axes,[[1,2,np.nan,4,5],[-.1,.1,np.nan,.2,.3]]):
            np.testing.assert_array_equal(track.lines[0].get_xdata(),table.pos)
            np.testing.assert_allclose(track.lines[0].get_ydata(),expected,equal_nan=True)
            self.assertEqual(len(track.collections[0].get_paths()),2)

    def test_no_flanks_features_cover_inclusive_gene_body(self):
        cfg,_ = coordinate_inputs(self.root)
        cfg['plot'].update(upstream=0,downstream=0)
        obj = FancyGene(cfg)
        obj.plot()
        rectangles = [p for p in obj.ax.patches if isinstance(p,Rectangle)][:3]
        self.assertAlmostEqual(rectangles[0].get_x(),.1)
        self.assertAlmostEqual(rectangles[-1].get_x()+rectangles[-1].get_width(),.9)
        self.assertEqual(obj.ax.child_axes[0].get_xlim(),(99.5,200.5))

    def test_cli_flanks_and_single_base_from_another_directory(self):
        other = self.root/'other'
        other.mkdir()
        env = dict(os.environ)
        if os.environ.get('HASTAT_CLI'):
            cli = [os.environ['HASTAT_CLI']]
        else:
            cli = [sys.executable,'-m','hastat.hastat']
            env['PYTHONPATH'] = str(Path(__file__).resolve().parents[1])
        for case in ['plus_flanks','minus_flanks','single_base']:
            cfg,_ = coordinate_inputs(self.root/case,case)
            path = self.root/case/'plot.pdf'
            result = subprocess.run(cli+['plot','gene','--gff',cfg['data']['gff_file'],'--genes','G1',
                '--toml',cfg['data']['toml_file'],'--upstream',str(cfg['plot']['upstream']),
                '--downstream',str(cfg['plot']['downstream']),'-o',str(path)],env=env,cwd=other,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertGreater(path.stat().st_size,1000)


if __name__ == '__main__':
    unittest.main(verbosity=2)
