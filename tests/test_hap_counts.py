"""Sample counts on gene haplotype rows, including real CLI integration."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from hastat.viz.gene2 import FancyGene
from test_gene_tracks import gene_gff
from test_population import write_vcf


def example_inputs(root, counts=(25, 12, 3)):
    root.mkdir(parents=True, exist_ok=True)
    names = ['001', '1', 'NA'] + [f's{i}' for i in range(4, sum(counts)+1)]
    labels = [f'Hap{i+1}' for i, count in enumerate(counts) for _ in range(count)]
    pd.DataFrame({'samples': names, 'haplotypes': labels}).to_csv(root / 'gene.group.csv', index=False)
    hap = pd.DataFrame({'chrom': ['1']*4, 'pos': [5, 10, 20, 25], 'ref': ['A']*4, 'alt': ['G']*4,
                        'Hap1': [0, 0, 2, 2], 'Hap2': [0, 2, 2, 0], 'Hap3': [2, 2, 0, 0]})
    hap.to_csv(root / 'gene.table.csv', index=False)
    gene_gff(root / 'gene.gff3')
    (root / 'hap.categories.csv').write_text('Hap1,A\nHap2,B\nHap3,B\n')
    return hap


def example_config(root, extra='', counts=True):
    if not counts:
        (root / 'gene.group.csv').unlink(missing_ok=True)
    (root / 'gene.toml').write_text('[default]\nbbox=[0.1,0.1,0.8,0.8]\nxy=[0.0,0.85]\n'
        'height=0.1\nhap_y=0.0\nhap_height=0.4\nhap_freq_show=false\n'
        'hap_file="gene"\n' + extra)
    return {'data': {'gff_file': str(root / 'gene.gff3'), 'genes': ['G1'], 'toml_file': str(root / 'gene.toml')},
            'plot': {'width': 8, 'height': 5}}


def labels_for(obj):
    return [text.get_text() for text in obj.ax.texts if text.get_text().startswith('Hap')]


class HaplotypeCountValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.addCleanup(plt.close, 'all')
        self.root = Path(self.tmp.name)
        example_inputs(self.root, (5, 3, 1))

    def draw(self, extra='', counts=True):
        obj = FancyGene(example_config(self.root, extra, counts))
        obj.plot()
        return obj

    def cli(self, *args, cwd=None):
        command = ([os.environ['HASTAT_CLI']] if os.environ.get('HASTAT_CLI')
                   else [sys.executable, '-m', 'hastat.hastat'])
        env = dict(os.environ)
        if not os.environ.get('HASTAT_CLI'):
            env['PYTHONPATH'] = str(Path(__file__).resolve().parents[1])
        return subprocess.run(command + list(map(str, args)), env=env, cwd=cwd, capture_output=True, text=True)

    def test_exact_sample_counts_with_literal_ids(self):
        self.assertEqual(labels_for(self.draw()), ['Hap1 (n=5)', 'Hap2 (n=3)', 'Hap3 (n=1)'])

    def test_no_mapping_keeps_existing_labels(self):
        self.assertEqual(labels_for(self.draw(counts=False)), ['Hap1','Hap2','Hap3'])

    def test_count_toggle_and_hidden_labels(self):
        self.assertEqual(labels_for(self.draw('hap_count_show=false\n')), ['Hap1','Hap2','Hap3'])
        self.assertEqual(labels_for(self.draw('hap_label_show=false\n')), [])
        with self.assertRaisesRegex(ValueError, 'true or false'):
            self.draw('hap_count_show="true"\n')

    def test_missing_group_is_optional_even_when_counts_enabled(self):
        self.assertEqual(labels_for(self.draw('hap_count_show=true\n', counts=False)), ['Hap1','Hap2','Hap3'])

    def test_duplicate_sample_ids_rejected_even_when_haplotype_hidden(self):
        path = self.root / 'gene.group.csv'
        path.write_text(path.read_text() + 's9,Hap3\n')
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            self.draw('hap_cols=["Hap1"]\n')

    def test_missing_displayed_assignment_is_not_invented_zero(self):
        path = self.root / 'gene.group.csv'
        path.write_text('samples,haplotypes\n001,Hap1\n1,Hap2\n')
        with self.assertRaisesRegex(ValueError, 'Hap3'):
            self.draw()
        with self.assertRaisesRegex(ValueError, 'Hap3'):
            self.draw('hap_count_show=false\n')
        self.assertEqual(labels_for(self.draw('hap_cols=["Hap2","Hap1"]\n')),
                         ['Hap1 (n=1)', 'Hap2 (n=1)'])

    def test_display_subset_keeps_original_counts(self):
        self.assertEqual(labels_for(self.draw('hap_cols=["Hap3","Hap1"]\n')),
                         ['Hap1 (n=5)', 'Hap3 (n=1)'])

    def test_visual_filter_does_not_recount_or_merge_haplotypes(self):
        path = self.root / 'gene.table.csv'
        table = pd.read_csv(path)
        table['Hap2'] = 1
        table.to_csv(path,index=False)
        self.assertEqual(labels_for(self.draw('hap_het_rate=[0.0,0.0]\n')),
                         ['Hap1 (n=5)','Hap3 (n=1)'])

    def test_group_sorting_retains_name_count_mapping(self):
        (self.root / 'hap.categories.csv').write_text('Hap1,Z\nHap2,A\nHap3,A\n')
        self.assertEqual(labels_for(self.draw('hap_color_file="hap.categories.csv"\n')),
                         ['Hap2 (n=3)','Hap3 (n=1)','Hap1 (n=5)'])

    def test_cluster_labels_and_bars_do_not_overlap(self):
        obj = self.draw('hap_color_file="hap.categories.csv"\nhap_cluster_show=true\n')
        self.assertEqual(set(labels_for(obj)), {'Hap1 (n=5)','Hap2 (n=3)','Hap3 (n=1)'})
        obj.fig.canvas.draw()
        renderer = obj.fig.canvas.get_renderer()
        texts, bars, tree = obj._count_label_layout[0]
        self.assertTrue(bars and tree)
        left = min(text.get_window_extent(renderer).x0 for text in texts)
        right = max(bar.get_window_extent(renderer).x1 for bar in bars)
        self.assertLess(right, left)
        self.assertLess(max(line.get_window_extent(renderer).x1 for line in tree),
                        min(bar.get_window_extent(renderer).x0 for bar in bars))
        # Export at a different DPI: extents must remain separated.
        obj.fig.set_dpi(200)
        obj.fig.canvas.draw()
        renderer = obj.fig.canvas.get_renderer()
        self.assertLess(max(bar.get_window_extent(renderer).x1 for bar in bars),
                        min(text.get_window_extent(renderer).x0 for text in texts))

    def test_cluster_without_groups_and_repeated_plot(self):
        obj = self.draw('hap_cluster_show=true\nhap_cluster_linkage="average"\n')
        self.assertEqual(set(labels_for(obj)), {'Hap1 (n=5)','Hap2 (n=3)','Hap3 (n=1)'})
        texts, bars, tree = obj._count_label_layout[0]
        obj.fig.canvas.draw()
        renderer = obj.fig.canvas.get_renderer()
        self.assertLess(max(line.get_window_extent(renderer).x1 for line in tree),
                        min(text.get_window_extent(renderer).x0 for text in texts))
        obj.plot()
        self.assertEqual(len(obj._count_label_layout), 1)
        self.assertEqual(len(labels_for(obj)), 3)

    def test_real_view_table_group_to_gene_plot_cli(self):
        samples = ('001','1','NA','s4','s5','s6')
        rows = [(5,'A','G',['0/0']*3+['1/1']*3),
                (15,'A','G',['0/0']*3+['1/1']*2+['0/0']),
                (25,'A','G',['1/1']*3+['0/0']*3)]
        vcf = write_vcf(self.root/'input.vcf',rows,samples=samples,chrom='1')
        for kind in ('table','group'):
            result = self.cli('view','-v',vcf,'-r','1:1-30','-t',kind,'-o',self.root/'gene')
            self.assertEqual(result.returncode,0,result.stderr)
        config = example_config(self.root)
        obj = FancyGene(config)
        obj.plot()
        self.assertEqual(set(labels_for(obj)), {'Hap1 (n=3)','Hap2 (n=2)','Hap3 (n=1)'})
        # Real console path, from outside the directory containing the inputs.
        other = self.root/'other'
        other.mkdir()
        output = self.root/'hap.svg'
        result = self.cli('plot','gene','--gff',self.root/'gene.gff3','--genes','G1',
                          '--toml',self.root/'gene.toml','-o',output,cwd=other)
        self.assertEqual(result.returncode,0,result.stderr)
        text = output.read_text()
        self.assertIn('Hap1 (n=3)',text)
        self.assertIn('Hap2 (n=2)',text)
        self.assertIn('Hap3 (n=1)',text)
        output.unlink()
        path = self.root/'gene.group.csv'
        path.write_text(path.read_text()+'001,Hap1\n')
        result = self.cli('plot','gene','--gff',self.root/'gene.gff3','--genes','G1',
                          '--toml',self.root/'gene.toml','-o',output,cwd=other)
        self.assertNotEqual(result.returncode,0)
        self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
