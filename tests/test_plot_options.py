"""Phenotype labels, omnibus tests, sample-size filtering and gene structures."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
import numpy as np
import pandas as pd
from scipy import stats
from hastat.viz.box import HapBox
from hastat.viz.gene import FancyGene
from hastat.utils.association import HaplotypeAssociation
from test_plot import frames, box_for
from test_hap_counts import example_inputs, example_config, labels_for
from test_population import write_vcf


class BoxOptions(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.addCleanup(plt.close, 'all')
        self.root = Path(tmp.name)

    def draw(self, method='WelchBonferroni', rotation=None):
        hap, phe = frames(([1, 2, 3, 4, np.nan], [3, 4, 5, 6, 7], [2, 4, 8, 16, 32]))
        hap.to_csv(self.root/'hap.csv', index=False)
        phe.to_csv(self.root/'phe.csv', index=False)
        config = {'data': {'sample_hap': str(self.root/'hap.csv'),
                           'sample_phe': str(self.root/'phe.csv')},
                  'plot': {'width': 6, 'height': 4, 'method': method, 'phe_index': 1,
                           'min_hap_size': 3, 'comparisons': [],
                           'save_fig': str(self.root/'box.png')}}
        if rotation is not None:
            config['plot']['x_rotation'] = rotation
        obj = HapBox(config)
        fig, ax = plt.subplots()
        obj.plot(ax=ax)
        return obj, ax, hap, phe

    def test_counts_use_current_trait_and_overall_without_brackets(self):
        obj, ax, hap, phe = self.draw()
        self.assertEqual([t.get_text() for t in ax.get_xticklabels()],
                         ['Hap1\n(n=4)', 'Hap2\n(n=5)', 'Hap3\n(n=5)'])
        self.assertTrue(all(t.get_rotation() == 0 for t in ax.get_xticklabels()))
        ref = HaplotypeAssociation(phe, hap, min_hap_size=3, method='WelchBonferroni')
        self.assertAlmostEqual(obj.omnibus_test('trait')['pvalue'], ref.anova_data['omnibus'][0]['pvalue'])
        self.assertIn('Welch ANOVA', ax.get_title())

    def test_rotation_keeps_counts_and_overall_annotation(self):
        for angle in (0, 45, 90, -37.5):
            with self.subTest(angle=angle):
                _, ax, _, _ = self.draw(rotation=angle)
                self.assertTrue(all(t.get_rotation() == angle % 360 for t in ax.get_xticklabels()))
                self.assertEqual(ax.get_xticklabels()[0].get_text(), 'Hap1\n(n=4)')
                self.assertIn('Welch ANOVA', ax.get_title())

    def test_overall_and_pairwise_share_full_input_family(self):
        obj = box_for(method='WelchBonferroni', labels=['Hap1', 'Hap2'])
        full = obj.omnibus_test('trait')
        self.assertEqual(full['groups'], 3)
        obj.add_comparisons([('Hap1', 'Hap2')])
        self.assertEqual(obj.omnibus_test('trait'), full)
        hap, phe = frames()
        ref = HaplotypeAssociation(phe, hap, min_hap_size=2, method='WelchBonferroni')
        self.assertAlmostEqual(full['pvalue'], ref.anova_data['omnibus'][0]['pvalue'])

    def test_ordinary_anova_matches_independent_scipy(self):
        arrays = ([1, 2, 3, 4, 5], [3, 4, 5, 6, 7], [1, 4, 7, 10, 13])
        for method in ('t-test', 'TukeyHSD'):
            obj = box_for(arrays, method=method)
            self.assertAlmostEqual(obj.omnibus_test('trait')['pvalue'], stats.f_oneway(*arrays).pvalue)

    def test_kruskal_known_reference_and_tie_correction(self):
        obj = box_for(([1, 3, 5, 7, 9], [2, 4, 6, 8, 10]), method='mannwhitneyu')
        self.assertAlmostEqual(obj.omnibus_test('trait')['pvalue'], .6015081344405895)
        tied = box_for(([1, 1, 1, 2, 2], [2, 2, 3, 3, 3]), method='mannwhitneyu')
        # Hand-calculated pooled ranks give H=5.4 after the tie correction.
        self.assertAlmostEqual(tied.omnibus_test('trait')['pvalue'], stats.chi2.sf(5.4, 1))
        small = box_for(([1, 2, 3], [4, 5, 6]), method='mannwhitneyu')
        self.assertIsNone(small.omnibus_test('trait'))

    def test_undefined_overall_and_extreme_p_formatting(self):
        self.assertIsNone(box_for(([1]*5, [2]*5), method='welch').omnibus_test('trait'))
        self.assertIsNone(box_for(([1, 2], [3, 4]), minimum=3).omnibus_test('trait'))
        self.assertIn('10^{-324}', HapBox._p_label(np.nextafter(0., 1.)))
        self.assertIn('P <', HapBox._p_label(0.))


class GeneOptions(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.addCleanup(plt.close, 'all')
        self.root = Path(tmp.name)

    def test_size_filter_changes_rows_but_preserves_af_and_ld(self):
        table = example_inputs(self.root, (5, 3, 1))
        samples = pd.read_csv(self.root/'gene.group.csv', dtype=str, keep_default_na=False)
        gt = {0:'0/0', 1:'0/1', 2:'1/1'}
        rows = [(r.pos, r.ref, r.alt, [gt[getattr(r, h)] for h in samples.haplotypes])
                for r in table.itertuples(index=False)]
        vcf = write_vcf(self.root/'test.vcf', rows, samples=samples.samples.tolist(), chrom='1')
        settings = f'vcf_file="{Path(vcf).as_posix()}"\nld_show=true\n'
        base_cfg = example_config(self.root, settings)
        p = self.root/'gene.toml'
        p.write_text(p.read_text().replace('af_show=false', 'af_show=true'))
        obj = FancyGene(base_cfg)
        obj.add_data(base_cfg['data'])
        conf = obj._get_gene_config('G1')
        gene = obj._extract_gene_features('G1')
        original = obj._prepare_haplotypes(gene, conf)
        filtered = obj._prepare_haplotypes(gene, dict(conf, hap_min_size=3, hap_count_show=False))
        self.assertEqual(filtered['haps'], ['Hap1', 'Hap2'])
        self.assertEqual(filtered['labels'], {'Hap1':'Hap1', 'Hap2':'Hap2'})
        pd.testing.assert_frame_equal(original['frequencies'], filtered['frequencies'])
        np.testing.assert_allclose(original['ld'], filtered['ld'], equal_nan=True)
        pd.testing.assert_frame_equal(original['sites'], filtered['sites'])

    def test_zero_filter_missing_counts_and_all_rows_hidden(self):
        example_inputs(self.root, (5, 3, 1))
        obj = FancyGene(example_config(self.root, 'hap_min_size=0\n'))
        obj.plot()
        self.assertEqual(len(labels_for(obj)), 3)
        obj = FancyGene(example_config(self.root, 'hap_min_size=99\n'))
        obj.plot()
        self.assertEqual(labels_for(obj), [])
        obj = FancyGene(example_config(self.root, 'hap_min_size=3\n', counts=False))
        with patch('hastat.viz.gene.logger.warning') as warning:
            obj.plot()
        self.assertTrue(warning.called)
        self.assertEqual(len(labels_for(obj)), 3)

    def gene(self, body, strand='+', settings=''):
        example_inputs(self.root)
        cfg = example_config(self.root, settings)
        (self.root/'gene.gff3').write_text(
            f'1\ttest\tgene\t1\t30\t.\t{strand}\t.\tID=G1\n' + body)
        obj = FancyGene(cfg)
        obj.plot()
        return obj

    def test_auto_cds_priority_and_unknown_feature_setting(self):
        body = ('1\ttest\texon\t1\t30\t.\t+\t.\tID=E1;Parent=G1\n'
                '1\ttest\tCDS\t5\t25\t.\t+\t0\tID=C1;Parent=G1\n')
        with patch('hastat.viz.gene.logger.warning') as warning:
            obj = self.gene(body, settings='gene_feature="exon"\n')
        self.assertTrue(any(call.args[1:] == ('gene_feature', 'default') for call in warning.call_args_list))
        self.assertNotIn('gene_feature', obj._get_gene_config('G1'))
        blocks = [p.get_gid() or '' for p in obj.ax.patches]
        self.assertEqual(sum(k.startswith('gene-CDS') for k in blocks), 1)
        self.assertFalse(any(k.startswith('gene-exon') for k in blocks))
        self.assertEqual(sum(k.startswith('gene-UTR') for k in blocks), 2)

    def test_gene_only_annotation_keeps_interval_and_direction(self):
        with patch('hastat.viz.gene.logger.warning') as warning:
            obj = self.gene('', strand='-')
        self.assertTrue(warning.called)
        self.assertEqual(sum(isinstance(p, PathPatch) for p in obj.ax.patches), 1)
        self.assertFalse(any((p.get_gid() or '').startswith('gene-') for p in obj.ax.patches))
        marker = next(t for t in obj.ax.texts if t.get_gid() == 'strand-direction')
        self.assertEqual(marker.get_text(), '<')

    def test_inferred_utrs_and_exon_only_fallback(self):
        body = ('1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=T1;Parent=G1\n'
                '1\ttest\texon\t1\t10\t.\t+\t.\tID=E1;Parent=T1\n'
                '1\ttest\texon\t21\t30\t.\t+\t.\tID=E2;Parent=T1\n')
        obj = self.gene(body + '1\ttest\tCDS\t5\t10\t.\t+\t0\tID=C1;Parent=T1\n'
                        '1\ttest\tCDS\t21\t25\t.\t+\t0\tID=C2;Parent=T1\n')
        self.assertEqual(obj.gene_data['UTR'], [[1, 4], [26, 30]])
        blocks = {p.get_gid():p for p in obj.ax.patches if p.get_gid()}
        cds = [p for k,p in blocks.items() if k.startswith('gene-CDS')]
        utr = [p for k,p in blocks.items() if k.startswith('gene-UTR')]
        self.assertEqual(len(utr), 2)
        self.assertAlmostEqual(utr[0].get_height()/cds[0].get_height(), .55)
        (self.root/'gene.gff3.sqlite3').unlink()
        obj = self.gene(body)
        self.assertEqual(obj.gene_data['UTR'], [])
        self.assertEqual(sum((p.get_gid() or '').startswith('gene-exon') for p in obj.ax.patches), 2)

    def test_generic_utr_and_transcript_separation(self):
        body = ('1\ttest\tmRNA\t1\t10\t.\t+\t.\tID=T1;Parent=G1\n'
                '1\ttest\texon\t1\t10\t.\t+\t.\tID=E1;Parent=T1\n'
                '1\ttest\tmRNA\t21\t30\t.\t+\t.\tID=T2;Parent=G1\n'
                '1\ttest\texon\t21\t30\t.\t+\t.\tID=E2;Parent=T2\n'
                '1\ttest\tCDS\t21\t25\t.\t+\t0\tID=C1;Parent=T2\n'
                '1\ttest\tUTR\t27\t30\t.\t+\t.\tID=U1;Parent=T2\n')
        obj = self.gene(body)
        # The exon-only isoform T1 must not acquire a fabricated UTR.
        self.assertEqual(obj.gene_data['UTR'], [[26, 30]])

    def test_partial_cds_does_not_infer_utr(self):
        body = ('1\ttest\tmRNA\t1\t30\t.\t-\t.\tID=T1;Parent=G1;partial=true\n'
                '1\ttest\texon\t1\t30\t.\t-\t.\tID=E1;Parent=T1\n'
                '1\ttest\tCDS\t5\t25\t.\t-\t0\tID=C1;Parent=T1\n')
        obj = self.gene(body, strand='-')
        self.assertEqual(obj.gene_data['UTR'], [])
        marker = next(t for t in obj.ax.texts if t.get_gid() == 'strand-direction')
        self.assertEqual(marker.get_text(), '<')
        self.assertEqual(marker.get_position(), (-.5, 0))


if __name__ == '__main__':
    unittest.main(verbosity=2)
