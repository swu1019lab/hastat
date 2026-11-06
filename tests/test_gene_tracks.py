"""Gene tracks: exact artist data and real view-to-plot CLI tests, PBS only."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import warnings
from xml.etree import ElementTree

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from hastat.dataio.table import read_window_table
from hastat.viz.gene import FancyGene
from test_population import write_vcf


def gene_gff(path, chrom='1'):
    path.write_text(f'{chrom}\ttest\tgene\t1\t30\t.\t+\t.\tID=G1\n'
                    f'{chrom}\ttest\tCDS\t3\t25\t.\t+\t0\tID=C1;Parent=G1\n')
    return path


class GeneTrackValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.addCleanup(plt.close, 'all')
        self.root = Path(self.tmp.name)
        self.gene = dict(chrom='1', start=1, end=30)

    def write(self, name, text):
        path = self.root / name
        path.write_text(text, encoding='utf-8')
        return path

    def obj(self):
        obj = FancyGene({'plot': {}})
        obj.fig, obj.ax = plt.subplots()
        return obj

    def cli(self, *args, cwd=None):
        command = ([os.environ['HASTAT_CLI']] if os.environ.get('HASTAT_CLI')
                   else [sys.executable, '-m', 'hastat.hastat'])
        env = dict(os.environ)
        if not os.environ.get('HASTAT_CLI'):
            env['PYTHONPATH'] = str(Path(__file__).resolve().parents[1])
        return subprocess.run(command + list(map(str, args)), cwd=cwd, env=env, capture_output=True, text=True)

    def test_literal_chromosome_and_group_labels(self):
        path = self.write('pi.csv', 'chrom,start,end,pi,groups\n01,1,5,0.1,NA\n1,1,5,0.2,01\nNA,1,5,,P\n')
        data = read_window_table(path, 'pi')
        self.assertEqual(data.chrom.tolist(), ['01', '1', 'NA'])
        self.assertEqual(data.groups.tolist(), ['NA', '01', 'P'])
        self.assertTrue(np.isnan(data.pi.iloc[2]))

    def test_invalid_schema_coordinates_and_duplicate_windows(self):
        texts = ['start,end,pi,pi\n1,2,0,0\n', 'start,end\n1,2\n',
                 'start,end,pi\n0,2,0\n', 'start,end,pi\n2,1,0\n',
                 'start,end,pi\n1.5,2,0\n', 'start,end,pi\n,2,0\n',
                 'start,end,pi\n1,inf,0\n', 'chrom,start,end,pi\n,1,2,0\n',
                 'start,end,pi,groups\n1,2,0,\n', 'start,end,pi\n1,2,0\n1,2,0.1\n']
        for text in texts:
            with self.subTest(text=text), self.assertRaises(ValueError):
                read_window_table(self.write('bad.csv', text), 'pi')

    def test_overlaps_and_same_window_in_distinct_groups_allowed(self):
        path = self.write('pi.csv', 'chrom,start,end,pi,groups\n1,1,10,.1,P\n1,5,15,.2,P\n1,1,10,.3,Q\n')
        self.assertEqual(len(read_window_table(path, 'pi')), 3)

    def test_invalid_values_rejected_and_undefined_retained(self):
        for metric, value in (('pi', '-.1'), ('pi', 'bad'), ('fst', 'inf'), ('fst', '-Infinity')):
            with self.subTest(metric=metric, value=value), self.assertRaises(ValueError):
                read_window_table(self.write('bad.csv', f'start,end,{metric}\n1,5,{value}\n'), metric)
        frame = read_window_table(self.write('fst.csv', 'start,end,fst\n1,5,-.2\n6,10,NA\n11,15,NaN\n'), 'fst')
        np.testing.assert_allclose(frame.fst, [-.2, np.nan, np.nan], equal_nan=True)

    def test_track_artist_values_and_order_preserve_na_group(self):
        path = self.write('pi.csv', 'chrom,start,end,pi,groups\n1,11,20,.002,NA\n1,1,10,.001,NA\n01,1,10,.9,P\n')
        obj = self.obj()
        lines = obj._draw_pi_plot(self.gene, {'pi_file': path}, [0, 0, 1, 1])
        self.assertEqual(len(lines), 1)
        self.assertEqual(lines[0].get_label(), 'NA')
        np.testing.assert_array_equal(lines[0].get_xdata(), [5.5, 15.5])
        np.testing.assert_allclose(lines[0].get_ydata(), [1, 2])
        self.assertEqual(obj.ax.child_axes[0].get_xlim(), (.5, 20.5))

    def test_single_window_has_marker_and_nonsingular_axis(self):
        obj = self.obj()
        path = self.write('pi.csv', 'chrom,start,end,pi\n1,5,10,.002\n')
        with warnings.catch_warnings(record=True) as caught:
            lines = obj._draw_pi_plot(self.gene, {'pi_file': path}, [0, 0, 1, 1])
        self.assertFalse(any('identical' in str(item.message).lower() for item in caught))
        self.assertEqual(lines[0].get_marker(), 'o')
        np.testing.assert_allclose(lines[0].get_xdata(), [7.5])
        np.testing.assert_allclose(lines[0].get_ydata(), [2])

    def test_single_base_track_axis(self):
        obj = self.obj()
        path = self.write('pi.csv', 'start,end,pi\n5,5,.1\n')
        obj._draw_pi_plot(dict(chrom='1', start=5, end=5), {'pi_file': path}, [0, 0, 1, 1])
        self.assertEqual(obj.ax.child_axes[0].get_xlim(), (4.5, 5.5))

    def test_nan_values_break_lines_and_fills(self):
        obj = self.obj()
        path = self.write('fst.csv', 'start,end,fst\n1,5,.1\n6,10,.2\n11,15,\n16,20,.3\n21,25,.4\n')
        lines = obj._draw_fst_plot(self.gene, {'fst_file': path}, [0, 0, 1, 1])
        np.testing.assert_allclose(lines[0].get_ydata(), [.1, .2, np.nan, .3, .4], equal_nan=True)
        self.assertEqual(len(obj.ax.child_axes[0].collections[0].get_paths()), 2)

    def test_all_undefined_not_drawn_as_zero(self):
        obj = self.obj()
        path = self.write('fst.csv', 'start,end,fst\n1,5,NaN\n6,10,\n')
        with self.assertLogs('hastat', level='WARNING'):
            lines = obj._draw_fst_plot(self.gene, {'fst_file': path}, [0, 0, 1, 1])
        self.assertTrue(np.isnan(lines[0].get_ydata()).all())
        self.assertEqual(len(obj.ax.child_axes[0].collections[0].get_paths()), 0)

    def test_negative_fst_and_line_style(self):
        obj = self.obj()
        path = self.write('fst.csv', 'chrom,start,end,fst,groups\n1,1,10,-.25,01_NA\n1,11,20,.1,01_NA\n')
        lines = obj._draw_fst_plot(self.gene, {'fst_file': path, 'fst_style': 'line'}, [0, 0, 1, 1])
        np.testing.assert_allclose(lines[0].get_ydata(), [-.25, .1])
        self.assertLess(obj.ax.child_axes[0].get_ylim()[0], -.25)
        self.assertEqual(len(obj.ax.child_axes[0].collections), 0)

    def test_window_track_uses_all_windows_on_target_chromosome(self):
        obj = self.obj()
        path = self.write('pi.csv', 'chrom,start,end,pi\n1,1,10,.1\n1,1,31,.2\n'
                          '1,100,150,\n2,1,500,.3\n')
        lines = obj._draw_pi_plot(self.gene, {'pi_file': path}, [0, 0, 1, 1])
        np.testing.assert_allclose(lines[0].get_xdata(), [5.5, 16, 125])
        np.testing.assert_allclose(lines[0].get_ydata(), [100, 200, np.nan], equal_nan=True)
        self.assertEqual(lines[0].axes.get_xlim(), (.5, 150.5))
        with self.assertLogs('hastat', level='WARNING'):
            self.assertIsNone(obj._draw_pi_plot(dict(self.gene, chrom='3'), {'pi_file': path}, [0,0,1,1]))

    def test_window_tracks_without_chrom_use_own_extent_and_clip_gene_shading(self):
        obj = self.obj()
        gene = dict(self.gene, gene_start=9, gene_end=20)
        for stat, bounds, shade in [('pi', (1, 11), (8.5, 3)), ('fst', (4, 42), (8.5, 12)),
                                    ('tajima_d', (40, 99), None)]:
            with self.subTest(stat=stat):
                lo, hi = bounds
                path = self.write(f'{stat}.csv', f'start,end,{stat}\n{lo},{hi},0.1\n')
                lines = obj._draw_window_track(gene, {stat+'_file': path, 'shade_gene': True}, [0,0,1,1], stat)
                axis = lines[0].axes
                self.assertEqual(axis.get_xlim(), (lo-.5, hi+.5))
                spans = [p for p in axis.patches if p.get_gid() == 'gene-span']
                self.assertEqual([(p.get_x(), p.get_width()) for p in spans], [shade] if shade else [])

    def test_empty_absent_and_invalid_style(self):
        obj = self.obj()
        self.assertIsNone(obj._draw_pi_plot(self.gene, {}, [0, 0, 1, 1]))
        path = self.write('empty.csv', 'start,end,pi\n')
        self.assertIsNone(obj._draw_pi_plot(self.gene, {'pi_file': path}, [0, 0, 1, 1]))
        path.write_text('start,end,pi\n1,5,.1\n')
        with self.assertRaisesRegex(ValueError, 'style'):
            obj._draw_pi_plot(self.gene, {'pi_file': path, 'pi_style': 'typo'}, [0, 0, 1, 1])

    def test_config_paths_resolve_relative_to_toml(self):
        path = self.write('pi.csv', 'start,end,pi\n1,5,.1\n')
        obj = self.obj()
        obj.toml_file = self.root / 'plot.toml'
        obj.toml_config = {'default': {'pi_file': 'pi.csv'}, 'G1': {'pi_style': 'line'}}
        config = obj._get_gene_config('G1')
        self.assertEqual(Path(config['pi_file']), path)
        self.assertEqual(config['pi_style'], 'line')
        self.assertEqual(obj.toml_config['default']['pi_file'], 'pi.csv')

    def test_missing_configured_files_raise_and_unreadable_inputs_raise(self):
        obj = self.obj()
        obj.toml_file = self.root / 'plot.toml'
        for key in ('pi_file', 'fst_file', 'tajima_d_file', 'ann_file', 'hap_file', 'hap_color_file'):
            obj.toml_config = {'default': {key: 'missing.csv'}}
            with self.subTest(key=key), self.assertRaises(FileNotFoundError):
                obj._get_gene_config('G1')
        bad = self.write('bad.csv', 'onecolumn\n')
        for method in (obj._load_annotation_data, obj._load_haplotype_group_data, obj._load_pi_data, obj._load_fst_data, obj._load_tajima_d_data):
            with self.subTest(method=method.__name__), self.assertRaises(ValueError):
                method(bad)
        with self.assertRaises(ValueError):
            obj._load_haplotype_data(bad, ['Hap1'])

    def test_numeric_chromosome_survives_annotation_and_haplotype_loading(self):
        obj = self.obj()
        ann = obj._load_annotation_data(self.write('ann.csv', 'chrom,pos,ann\n01,5,missense\n'))
        hap = obj._load_haplotype_data(self.write('hap.csv', 'chrom,pos,ref,alt,Hap1\n01,5,A,G,0\n'))
        self.assertEqual(ann.chrom.iloc[0], '01')
        self.assertEqual(hap.chrom.iloc[0], '01')

    def test_real_view_to_gene_plot_cli_from_another_directory(self):
        samples = ('001', '1', 'NA', 's4')
        rows = [(5, 'A', 'G', ['0/0', '0/1', '1/1', '0/1']),
                (15, 'A', 'G', ['1/1', '0/1', '0/0', '0/1']),
                (25, 'A', 'G', ['0/1', '0/1', '0/1', '0/1'])]
        vcf = write_vcf(self.root / 'input.vcf', rows, samples=samples, chrom='1')
        group = self.write('groups.csv', 'samples,groups\n001,NA\n1,NA\nNA,01\ns4,01\n')
        for metric in ('pi', 'fst', 'tajima_d'):
            result = self.cli('view', '-v', vcf, '-r', '1:1-30', '-g', group, '-t', metric,
                              '--size', 30, '--step', 30, '-o', self.root / 'windows')
            self.assertEqual(result.returncode, 0, result.stderr)
        gff = gene_gff(self.root / 'genes.gff3')
        config = self.write('gene.toml', '[default]\nbbox=[0.1,0.1,0.8,0.8]\npi_file="windows.pi.csv"\nfst_file="windows.fst.csv"\ntajima_d_file="windows.tajima_d.csv"\n')
        other = self.root / 'other'
        other.mkdir()
        output = self.root / 'gene.pdf'
        result = self.cli('plot', 'gene', '--gff', gff, '--gene', 'G1', '--toml', config, '-o', output, cwd=other)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertGreater(output.stat().st_size, 1000)
        self.assertEqual(set(read_window_table(self.root / 'windows.pi.csv', 'pi').groups), {'01', 'NA'})
        output.unlink()
        config.write_text('[default]\npi_file="missing.csv"\n')
        result = self.cli('plot', 'gene', '--gff', gff, '--gene', 'G1', '--toml', config, '-o', output, cwd=other)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(output.exists())

    def test_tajima_signed_values_nan_gap_and_zero_reference(self):
        obj = self.obj()
        path = self.write('d.csv', 'chrom,start,end,counts,tajima_d,groups\n1,1,5,3,-1.2,NA\n1,6,10,1,,NA\n1,11,15,3,2.1,NA\n')
        obj.toml_file = self.root / 'plot.toml'
        obj.toml_config = {'default': {'tajima_d_file': 'd.csv', 'tajima_d_style': 'line'}}
        config = obj._get_gene_config('G1')
        self.assertEqual(Path(config['tajima_d_file']), path)
        lines = obj._draw_tajima_d_plot(self.gene, config, [0, 0, 1, 1])
        np.testing.assert_allclose(lines[0].get_ydata(), [-1.2, np.nan, 2.1], equal_nan=True)
        np.testing.assert_allclose(lines[0].get_xdata(), [3, 8, 13])
        self.assertEqual(lines[0].get_label(), 'NA')
        axis = obj.ax.child_axes[0]
        self.assertEqual(axis.get_ylabel(), "Tajima's $D$")
        np.testing.assert_array_equal(axis.lines[0].get_ydata(), [0, 0])
        self.assertEqual(axis.get_xlim(), (.5, 15.5))
        self.assertIsNone(obj._draw_tajima_d_plot(self.gene, {}, [0, 0, 1, 1]))

    def test_tajima_all_undefined_keeps_nan(self):
        obj = self.obj()
        path = self.write('d.csv', 'start,end,tajima_d\n1,5,NaN\n6,10,\n')
        with self.assertLogs('hastat', level='WARNING'):
            lines = obj._draw_tajima_d_plot(self.gene, {'tajima_d_file': path}, [0, 0, 1, 1])
        self.assertTrue(np.isnan(lines[0].get_ydata()).all())
        self.assertEqual(len(obj.ax.child_axes[0].collections[0].get_paths()), 0)

    def test_gene_only_config_without_default_table(self):
        gff = gene_gff(self.root / 'genes.gff3')
        config = self.write('gene.toml', '[G1]\ngene_color="C1"\n')
        obj = FancyGene({'data': {'gff_file': str(gff), 'gene': 'G1', 'toml_file': str(config)}, 'plot': {}})
        fig, ax = obj.plot()
        self.assertGreater(len(ax.patches), 0)
        self.assertEqual(len(ax.child_axes), 0)

    def test_cli_selects_one_gene_per_file_and_ignores_unknown_keys(self):
        gff = gene_gff(self.root / 'two.gff3')
        with gff.open('a') as stream:
            stream.write('1\ttest\tgene\t51\t80\t.\t-\t.\tID=G2\n'
                         '1\ttest\texon\t54\t76\t.\t-\t.\tID=E2;Parent=G2\n')
        config = self.write('two.toml', '[default]\ngene_color="#154276"\n'
            '[G1]\nwidth=0.1\n[G2]\ngene_color="#ad534d"\ngene_feature="exon"\n')
        for gene, other, color in [('G1','G2','#154276'), ('G2','G1','#ad534d')]:
            output = self.root / f'{gene}.svg'
            result = self.cli('plot', 'gene', '--gff', gff, '--gene', gene, '--toml', config, '-o', output)
            self.assertEqual(result.returncode, 0, result.stderr)
            svg = ElementTree.parse(output)
            labels = [t.text for t in svg.iter('{http://www.w3.org/2000/svg}text')]
            self.assertIn(gene, labels)
            self.assertNotIn(other, labels)
            self.assertIn(f'fill: {color}', output.read_text())
            if gene == 'G1':
                self.assertIn("Unknown gene TOML parameter 'width' in [G1]", result.stderr)

    def test_cli_rejects_missing_plural_multiple_or_repeated_gene(self):
        output = self.root / 'invalid.svg'
        for selection in [[], ['--genes','G1'], ['--gene','G1','G2'], ['--gene','G1','--gene','G2']]:
            with self.subTest(selection=selection):
                result = self.cli('plot', 'gene', '--gff', 'unused.gff3', '--toml', 'unused.toml',
                                  *selection, '-o', output)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
