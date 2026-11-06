"""Sample joins, statistical annotations and real CLI regression; run inside PBS."""
from itertools import combinations
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
from scipy import stats

from hastat.viz.box import HapBox
from hastat.viz.bar import HapBar
from hastat.viz.pie import HapPie
from hastat.dataio.table import read_sample_csv, parse_numeric_values
from hastat.utils.association import HaplotypeAssociation

ROOT = Path(__file__).resolve().parents[1]
ARRAYS = ([1, 2, 3, 4, 5], [3, 4, 5, 6, 7], [1, 4, 7, 10, 13])


def frames(arrays=ARRAYS):
    rows = [(f's{i}_{j}', f'Hap{i+1}', value) for i, a in enumerate(arrays) for j, value in enumerate(a)]
    return (pd.DataFrame([(s, h) for s, h, _ in rows], columns=['samples', 'haplotypes']),
            pd.DataFrame([(s, y) for s, _, y in rows], columns=['samples', 'trait']))


def box_for(arrays=ARRAYS, method='t-test', labels=None, minimum=2):
    box = HapBox({'plot': {'min_hap_size': minimum}})
    hap, phe = frames(arrays)
    box.add_hap(hap, labels)
    box.add_phe(phe)
    box.add_method(method)
    return box


class PlotValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.addCleanup(plt.close, 'all')
        self.root = Path(self.tmp.name)

    def write(self, name, text):
        path = self.root / name
        path.write_text(text, encoding='utf-8')
        return str(path)

    def config(self, arrays=ARRAYS, method='t-test', comparisons=None):
        hap, phe = frames(arrays)
        hap.to_csv(self.root / 'hap.csv', index=False)
        phe.to_csv(self.root / 'phe.csv', index=False)
        return {'data': {'sample_hap': str(self.root / 'hap.csv'), 'sample_phe': str(self.root / 'phe.csv')},
                'plot': {'haplotypes': [], 'comparisons': comparisons or [('Hap1', 'Hap2')],
                         'method': method, 'phe_index': 1, 'step_size': [], 'width': 5, 'height': 4,
                         'save_fig': str(self.root / 'plot.pdf')}}

    def test_literal_ids_and_group_labels_survive_csv(self):
        hap = self.write('hap.csv', 'samples,haplotypes\n001,Hap1\n1,Hap2\nNA,Hap1\n')
        meta = self.write('meta.csv', 'samples,pop\n001,NA\n1,P\nNA,P\n')
        for cls in (HapBar, HapPie):
            obj = cls({'plot': {}})
            obj.add_data({'sample_hap': hap, 'sample_group': meta})
            result = obj.get_plot_data('pop')
            self.assertEqual(result.to_dict(), {'NA': {'Hap1': 1, 'Hap2': 0}, 'P': {'Hap1': 1, 'Hap2': 1}})
            self.assertEqual(set(obj.data['sample']), {'001', '1', 'NA'})
        phe = self.write('phe.csv', 'samples,trait\n001,0\n1,2\nNA,NA\n')
        obj = HapBox({'plot': {}})
        obj.add_data({'sample_hap': hap, 'sample_phe': phe})
        self.assertEqual(obj.get_plot_data('trait')[1].tolist(), [0, 2])
        self.assertEqual(set(obj.data['sample']), {'001', '1', 'NA'})

    def test_duplicate_hap_ids_rejected_before_filtering(self):
        frame = pd.DataFrame({'id': ['001', '001'], 'hap': ['A', 'B']})
        for cls in (HapBox, HapBar, HapPie):
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                cls({'plot': {}}).add_hap(frame, ['A'])

    def test_duplicate_metadata_ids_rejected(self):
        frame = pd.DataFrame({'id': ['001', ' 001 '], 'trait': [1, 2]})
        for cls in (HapBar, HapPie):
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                cls({'plot': {}}).add_group(frame)
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            HapBox({'plot': {}}).add_phe(frame)

    def test_missing_ids_labels_and_bad_shapes_rejected(self):
        for cls in (HapBox, HapBar, HapPie):
            for frame in (pd.DataFrame({'id': [''], 'hap': ['A']}),
                          pd.DataFrame({'id': ['x'], 'hap': [None]}),
                          pd.DataFrame({'id': ['x']}),
                          pd.DataFrame([['x', 'A', 'B']], columns=['id', 'a', 'b'])):
                with self.assertRaises(ValueError):
                    cls({'plot': {}}).add_hap(frame)

    def test_duplicate_csv_headers_and_reserved_names_rejected(self):
        with self.assertRaises(ValueError):
            read_sample_csv(self.write('bad.csv', 'id,trait,trait\nx,1,2\n'))
        obj = HapBox({'plot': {}})
        for columns in (['id', 'sample'], ['id', 'haplotypes']):
            with self.assertRaises(ValueError):
                obj.add_phe(pd.DataFrame([['x', 1]], columns=columns))

    def test_input_frames_unchanged(self):
        hap, phe = frames()
        before_hap, before_phe = hap.copy(), phe.copy()
        obj = HapBox({'plot': {}})
        obj.add_hap(hap)
        obj.add_phe(phe, 'trait')
        obj.get_plot_data('trait')
        pd.testing.assert_frame_equal(hap, before_hap)
        pd.testing.assert_frame_equal(phe, before_phe)

    def test_subset_join_and_counts_preserve_unique_samples(self):
        hap, _ = frames()
        group = pd.DataFrame({'id': ['s0_0', 's0_1', 's1_0', 'extra'], 'pop': ['P', 'Q', 'P', 'P']})
        for cls in (HapBar, HapPie):
            obj = cls({'plot': {}})
            obj.add_hap(hap)
            obj.add_group(group)
            counts = obj.get_plot_data('pop')
            self.assertEqual(int(counts.to_numpy().sum()), 3)
            np.testing.assert_allclose(obj.get_plot_data('pop', True).sum(), 1)

    def test_empty_or_nonoverlap_inputs_raise(self):
        for cls in (HapBox, HapBar, HapPie):
            obj = cls({'plot': {}})
            obj.add_hap(pd.DataFrame({'id': ['a'], 'hap': ['A']}))
            meta = pd.DataFrame({'id': ['b'], 'value': [1]})
            (obj.add_phe if cls is HapBox else obj.add_group)(meta)
            with self.assertRaisesRegex(ValueError, 'overlapping'):
                obj.get_plot_data('value')

    def test_numeric_missing_infinite_and_bad_text(self):
        values = parse_numeric_values(pd.Series(['0', '', 'NA', 'NaN', 'inf', '-inf', '-2']))
        self.assertEqual(values.dropna().tolist(), [0, -2])
        with self.assertRaises(ValueError):
            parse_numeric_values(pd.Series(['oops']))
        obj = box_for(([np.nan, np.inf], [np.nan, np.inf]))
        with self.assertRaisesRegex(ValueError, 'finite'):
            obj.get_plot_data('trait')

    def test_unestimable_never_means_ns(self):
        obj = box_for()
        for value in (np.nan, np.inf, -np.inf, None, -0.1, 1.1):
            self.assertIsNone(obj.get_significance_symbol(value))
        self.assertEqual([obj.get_significance_symbol(p) for p in (0, .001, .01, .05, 1)],
                         ['***', '**', '*', 'ns', 'ns'])
        self.assertTrue(np.isnan(obj.perform_statistical_test([1, 1], [2, 2])))
        self.assertTrue(np.isnan(obj.perform_statistical_test([1], [2, 3])))

    def test_pooled_t_and_full_family_bonferroni_formula(self):
        result = box_for().comparison_pvalues('trait')
        for i, j in combinations(range(3), 2):
            a, b = np.asarray(ARRAYS[i]), np.asarray(ARRAYS[j])
            df = len(a) + len(b) - 2
            pooled = ((len(a)-1)*a.var(ddof=1) + (len(b)-1)*b.var(ddof=1)) / df
            t = (a.mean() - b.mean()) / np.sqrt(pooled * (1/len(a) + 1/len(b)))
            reference = min(1, 3 * 2 * stats.t.sf(abs(t), df))
            self.assertAlmostEqual(result[f'Hap{i+1}', f'Hap{j+1}'], reference, places=13)

    def test_tukey_and_welch_match_stat_exactly(self):
        hap, phe = frames()
        for method in ('TukeyHSD', 'WelchBonferroni'):
            result = box_for(method=method).comparison_pvalues('trait')
            reference = HaplotypeAssociation(phe, hap, min_hap_size=2, method=method).comparisons
            for row in reference.to_dict('records'):
                self.assertEqual(result[row['group1'], row['group2']], row['p-adj'])
        self.assertEqual(box_for(method='welch').comparison_pvalues('trait'),
                         box_for(method='WelchBonferroni').comparison_pvalues('trait'))

    def test_tukey_independent_scipy_reference(self):
        result = box_for(method='TukeyHSD').comparison_pvalues('trait')
        reference = stats.tukey_hsd(*ARRAYS).pvalue
        for i, j in combinations(range(3), 2):
            self.assertAlmostEqual(result[f'Hap{i+1}', f'Hap{j+1}'], reference[i, j], places=12)

    def test_display_filters_and_brackets_do_not_change_pvalues(self):
        for method in ('t-test', 'mannwhitneyu', 'TukeyHSD', 'WelchBonferroni'):
            obj = box_for(method=method, labels=['Hap1', 'Hap2'])
            obj.add_comparisons([('Hap1', 'Hap2')])
            self.assertEqual(obj.comparison_pvalues('trait'), box_for(method=method).comparison_pvalues('trait'))

    def test_mann_whitney_exact_enumeration(self):
        a, b = [1, 3, 6], [2, 4, 5, 7]
        obj = box_for((a, b), method='mannwhitneyu')
        observed = sum(x > y for x in a for y in b)
        pooled = sorted(a+b)
        all_u = []
        for chosen in combinations(range(len(pooled)), len(a)):
            other = set(range(len(pooled))) - set(chosen)
            all_u.append(sum(pooled[i] > pooled[j] for i in chosen for j in other))
        ref = min(1, 2 * sum(u >= max(observed, len(a)*len(b)-observed) for u in all_u) / len(all_u))
        self.assertAlmostEqual(obj.comparison_pvalues('trait')['Hap1', 'Hap2'], ref, places=14)

    def test_mann_whitney_tie_corrected_normal_reference(self):
        a, b = [1, 1, 2, 3, 3], [2, 3, 4, 4, 5]
        u = sum((x > y) + 0.5*(x == y) for x in a for y in b)
        _, counts = np.unique(a+b, return_counts=True)
        n, m = len(a), len(b)
        variance = n*m/12 * (n+m+1 - sum(counts**3-counts)/((n+m)*(n+m-1)))
        z = (max(u, n*m-u) - n*m/2 - 0.5) / np.sqrt(variance)
        ref = min(1, 2*stats.norm.sf(z))
        self.assertAlmostEqual(box_for((a, b), method='mannwhitneyu').comparison_pvalues('trait')['Hap1', 'Hap2'], ref, places=14)

    def test_eligibility_and_zero_variance_match_stat(self):
        for method in ('TukeyHSD', 'WelchBonferroni'):
            self.assertEqual(box_for(([1], [2]), method=method).comparison_pvalues('trait'), {})
            self.assertEqual(box_for(([1, 1], [2, 2]), method=method).comparison_pvalues('trait'), {})
        self.assertEqual(len(box_for(minimum=6).comparison_pvalues('trait')), 0)

    def test_invalid_methods_pairs_and_minimum(self):
        obj = box_for()
        with self.assertRaises(ValueError):
            obj.add_method('AllPairTest')
        for pairs in ([('A', 'A')], [('A', 'B'), ('B', 'A')], [('A',)], ['AB']):
            with self.assertRaises(ValueError):
                obj.add_comparisons(pairs)
        for value in (1, 0, True, 2.5):
            with self.assertRaises(ValueError):
                HapBox({'plot': {'min_hap_size': value}})

    def test_plot_skips_invalid_annotation_and_keeps_points(self):
        obj = HapBox(self.config(([1, 1], [2, 2])))
        fig, ax = plt.subplots()
        obj.plot(ax)
        self.assertEqual(len(ax.texts), 0)
        self.assertEqual(len(ax.collections[0].get_offsets()), 4)
        self.assertTrue((self.root / 'plot.pdf').is_file())

    def test_negative_data_brackets_above_points_and_rng_unchanged(self):
        obj = HapBox(self.config(([-10, -9, -8], [-5, -4, -3])))
        state = np.random.get_state()
        fig, ax = plt.subplots()
        obj.plot(ax)
        self.assertGreater(ax.texts[0].get_position()[1], -3)
        current = np.random.get_state()
        np.testing.assert_array_equal(state[1], current[1])
        self.assertEqual(state[2:], current[2:])

    def test_real_cli_three_plot_types_and_failure_exit(self):
        cfg = self.config(method='TukeyHSD')
        hap, phe = frames()
        meta = hap.rename(columns={'haplotypes': 'pop'})
        meta['pop'] = ['P' if i % 2 else 'Q' for i in range(len(meta))]
        meta.to_csv(self.root / 'meta.csv', index=False)
        for kind in ('box', 'bar', 'pie'):
            command = [sys.executable, '-m', 'hastat.hastat', 'plot', kind,
                       '--sample_hap', cfg['data']['sample_hap'], '-o', str(self.root / f'{kind}.pdf')]
            command += (['--sample_phe', cfg['data']['sample_phe'], '--method', 'TukeyHSD',
                         '--comparisons', 'Hap1', 'Hap2'] if kind == 'box'
                        else ['--sample_group', str(self.root / 'meta.csv')])
            proc = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertGreater((self.root / f'{kind}.pdf').stat().st_size, 100)
            if kind == 'box':
                pd.concat([hap, hap.iloc[:1]]).to_csv(self.root / 'duplicate.csv', index=False)
                command[command.index('--sample_hap')+1] = str(self.root / 'duplicate.csv')
                command[command.index('-o')+1] = str(self.root / 'invalid.pdf')
                bad = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
                self.assertNotEqual(bad.returncode, 0)
                self.assertIn('duplicate', bad.stderr)
                self.assertFalse((self.root / 'invalid.pdf').exists())
        self.assertFalse(list(self.root.glob('*.summary.json')))


if __name__ == '__main__':
    unittest.main()
