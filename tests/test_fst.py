"""Independent WC84 and VCF/CLI integration checks; run inside PBS."""
from hastat.utils import population
import argparse
import ast
from importlib.util import find_spec
import itertools
import logging
import math
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import allel
import numpy as np
import pandas as pd
import pysam

from hastat.func import view
from hastat.utils.variant import VariantData
from test_population import write_vcf

ROOT = Path(__file__).resolve().parents[1]
SAMPLES = ('s1', 's2', 's3', 's4')
DIFFERENT = ['0/0', '0/0', '1/1', '1/1']
HETERO = ['0/1'] * 4


def wc_components(calls, group1, group2):
    """Scalar allele-1 WC84 components from individual calls, no library stats."""
    called = []
    for indices in (group1, group2):
        alleles = [str(calls[i]).replace('|', '/').split('/') for i in indices]
        called.append([tuple(map(int, pair)) for pair in alleles if len(pair) == 2 and '.' not in pair])
    n = [len(pop) for pop in called]
    if min(n) == 0 or sum(n) <= 2:
        return None
    total, mean_n = sum(n), sum(n) / 2
    nc = total - sum(x * x for x in n) / total
    p = [sum(sum(gt) for gt in pop) / (2 * len(pop)) for pop in called]
    pbar = sum(ni * pi for ni, pi in zip(n, p)) / total
    s2 = sum(ni * (pi - pbar) ** 2 for ni, pi in zip(n, p)) / mean_n
    h = sum(gt[0] != gt[1] for pop in called for gt in pop) / total
    a = mean_n / nc * (s2 - (pbar * (1 - pbar) - s2 / 2 - h / 4) / (mean_n - 1))
    b = mean_n / (mean_n - 1) * (pbar * (1 - pbar) - s2 / 2 - (2 * mean_n - 1) * h / (4 * mean_n))
    c = h / 2
    return a, b, c


def wc_window(rows, group1, group2):
    parts = [wc_components(gt, group1, group2) for _, _, _, gt in rows]
    parts = [part for part in parts if part is not None]
    if not parts:
        return math.nan
    a, b, c = [sum(part[i] for part in parts) for i in range(3)]
    return a / (a + b + c) if a + b + c != 0 else math.nan


class FstValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root, self.counter = Path(self.tmp.name), 0

    def variant(self, rows, samples=SAMPLES):
        self.counter += 1
        path = write_vcf(self.root / f'input{self.counter}.vcf', rows, samples=samples)
        handle = pysam.VariantFile(path)
        self.addCleanup(handle.close)
        return VariantData(handle), path

    def calculate(self, gv, start=1, end=20, size=20, step=20, groups=(SAMPLES[:2], SAMPLES[2:])):
        return population.calculate_fst(gv, 'chr1', start, end, list(groups[0]), list(groups[1]), size=size, step=step)

    def close(self, observed, expected):
        np.testing.assert_allclose(observed, expected, rtol=1e-11, atol=1e-12, equal_nan=True)

    def test_exhaustive_four_individual_states(self):
        states = ('0/0', '0/1', '1/1', './.')
        rows = [(i, 'A', 'G', list(gt)) for i, gt in enumerate(itertools.product(states, repeat=4), 1)]
        gv, _ = self.variant(rows)
        result = self.calculate(gv, 1, len(rows), 1, 1)
        self.assertEqual(len(result), 256)
        self.close(result[:, -1].astype(float), [wc_window([row], [0, 1], [2, 3]) for row in rows])
        self.assertTrue(np.all(result[:, -2].astype(int) == 1))

    def test_library_example_and_negative_value(self):
        rows = [(1, 'A', 'G', DIFFERENT), (2, 'A', 'G', HETERO),
                (3, 'A', 'G', ['0/0'] * 4), (4, 'A', 'G', ['0/0', '1/1', '0/1', './.'])]
        gv, _ = self.variant(rows)
        values = self.calculate(gv, 1, 4, 1, 1)[:, -1].astype(float)
        self.close(values, [1, 0, np.nan, -1.8])
        self.close(values[-1], wc_window(rows[-1:], [0, 1], [2, 3]))

    def test_ratio_of_sums_not_mean_ratios(self):
        rows = [(1, 'A', 'G', DIFFERENT), (2, 'A', 'G', HETERO)]
        gv, _ = self.variant(rows)
        value = float(self.calculate(gv, 1, 2, 2, 2)[0, -1])
        self.close(value, 2 / 3)
        self.assertNotAlmostEqual(value, .5)

    def test_seeded_unbalanced_missing_reference(self):
        rng = np.random.default_rng(20261031)
        names = tuple(f'u{i}' for i in range(9))
        states = np.array(['0/0', '0/1', '1/1', './.'])
        rows = [(i, 'A', 'G', rng.choice(states, 9, p=[.35, .25, .3, .1]).tolist()) for i in range(1, 301)]
        gv, _ = self.variant(rows, samples=names)
        result = self.calculate(gv, 1, 300, 11, 5, groups=(names[:3], names[3:]))
        for row in result:
            left, right = map(int, row[1:3])
            expected = wc_window([r for r in rows if left <= r[0] <= right], list(range(3)), list(range(3, 9)))
            self.close(float(row[-1]), expected)

    def test_unestimable_components_do_not_dilute_window(self):
        rows = [(1, 'A', 'G', DIFFERENT), (2, 'A', 'G', ['./.', './.', '0/1', '0/1']),
                (3, 'A', 'G', ['0/1', './.', '0/1', './.'])]
        gv, _ = self.variant(rows)
        result = self.calculate(gv, 1, 3, 3, 3)
        self.close(float(result[0, -1]), 1)
        self.assertEqual(int(result[0, -2]), 3)
        self.assertTrue(np.isnan(self.calculate(gv, 2, 3, 1, 1)[:, -1].astype(float)).all())

    def test_partial_calls_and_outside_population_missing(self):
        rows = [(1, 'A', 'G', ['0/.', '1/1', '0/0', '0/1', './.']),
                (2, 'A', 'G', ['./0', '1/1', '0/0', '0/1', '0/.'])]
        # VCF represents a partially missing reference call as ./0.
        gv, _ = self.variant(rows, samples=SAMPLES + ('outside',))
        values = self.calculate(gv, 1, 2, 1, 1)[:, -1].astype(float)
        self.close(values, [wc_window([r], [0, 1], [2, 3]) for r in rows])
        self.close(values[0], values[1])

    def test_one_plus_two_allowed_one_plus_one_undefined(self):
        rows = [(1, 'A', 'G', DIFFERENT)]
        gv, _ = self.variant(rows)
        self.close(float(self.calculate(gv, 1, 1, 1, 1, groups=(['s1'], ['s3', 's4']))[0, -1]), 1)
        self.assertTrue(math.isnan(float(self.calculate(gv, 1, 1, 1, 1, groups=(['s1'], ['s3']))[0, -1])))

    def test_empty_all_missing_and_monomorphic(self):
        for rows in ([], [(10, 'A', 'G', ['./.'] * 4)], [(10, 'A', 'G', ['0/0'] * 4)]):
            gv, _ = self.variant(rows)
            result = self.calculate(gv, 1, 20, 10, 10)
            self.assertEqual(result.shape, (2, 5))
            self.assertTrue(np.isnan(result[:, -1].astype(float)).all())
            self.assertEqual(result[:, -2].astype(int).sum(), len(rows))

    def test_tail_single_base_overlap_and_gaps(self):
        rows = [(p, 'A', 'G', DIFFERENT if p % 2 else HETERO) for p in range(1, 12)]
        gv, _ = self.variant(rows)
        for size, step, expected in ((5, 5, [[1, 5], [6, 10], [11, 11]]),
                                     (5, 3, [[1, 5], [4, 8], [7, 11]]),
                                     (2, 4, [[1, 2], [5, 6], [9, 10]])):
            result = self.calculate(gv, 1, 11, size, step)
            np.testing.assert_array_equal(result[:, 1:3].astype(int), expected)
            for output in result:
                left, right = map(int, output[1:3])
                self.close(float(output[-1]), wc_window([r for r in rows if left <= r[0] <= right], [0, 1], [2, 3]))
        self.close(float(self.calculate(gv, 11, 11, 1000, 1000)[0, -1]), 1)

    def test_invalid_population_and_window_parameters(self):
        gv, _ = self.variant([(1, 'A', 'G', DIFFERENT)])
        for groups in (([], ['s3']), (['absent'], ['s3']), (['s1', 's1'], ['s3']),
                       (['s1', 's2'], ['s2', 's3'])):
            with self.subTest(groups=groups), self.assertRaises(ValueError):
                self.calculate(gv, groups=groups)
        for samples in (None, 's1'):
            with self.assertRaises(ValueError):
                population.calculate_fst(gv, 'chr1', 1, 1, samples, ['s3'], size=1, step=1)
        for size, step in ((0, 1), (1, 0), (1.5, 1), (1, -1)):
            with self.assertRaises(ValueError):
                self.calculate(gv, size=size, step=step)

    def test_population_order_and_unused_samples(self):
        rows = [(1, 'A', 'G', DIFFERENT + ['./.']), (2, 'A', 'G', ['0/0', '1/1', '0/1', './.', '1/1'])]
        gv, _ = self.variant(rows, samples=SAMPLES + ('outside',))
        baseline = self.calculate(gv)[:, -1].astype(float)
        for groups in ((['s2', 's1'], ['s4', 's3']), (['s3', 's4'], ['s1', 's2'])):
            self.close(self.calculate(gv, groups=groups)[:, -1].astype(float), baseline)

    def test_numpy_settings_restored_on_success_and_error(self):
        gv, _ = self.variant([(1, 'A', 'G', ['0/0'] * 4)])
        original = np.geterr()
        try:
            np.seterr(all='raise')
            expected = np.geterr()
            self.calculate(gv)
            self.assertEqual(np.geterr(), expected)
            with patch('hastat.utils.population.allel.windowed_weir_cockerham_fst', side_effect=RuntimeError('injected')):
                with self.assertRaises(RuntimeError):
                    self.calculate(gv)
            self.assertEqual(np.geterr(), expected)
        finally:
            np.seterr(**original)

    def test_actual_library_called_with_explicit_windows_and_usable_sites(self):
        rows = [(1, 'A', 'G', DIFFERENT), (2, 'A', 'G', ['./.'] * 4), (11, 'A', 'G', HETERO)]
        gv, _ = self.variant(rows)
        real = allel.windowed_weir_cockerham_fst
        with patch('hastat.utils.population.allel.windowed_weir_cockerham_fst', wraps=real) as call:
            result = self.calculate(gv, 1, 11, 5, 5)
        self.assertEqual(call.call_count, 1)
        np.testing.assert_array_equal(call.call_args.args[0], [1, 11])
        np.testing.assert_array_equal(call.call_args.kwargs['windows'], [[1, 5], [6, 10], [11, 11]])
        self.close(result[:, -1].astype(float), [1, np.nan, 0])
        np.testing.assert_array_equal(result[:, -2].astype(int), [2, 0, 1])

    def test_biallelic_indel_scope_retained(self):
        gv, _ = self.variant([(1, 'A', 'AG', DIFFERENT)])
        result = self.calculate(gv)
        self.close(float(result[0, -1]), 1)
        self.assertEqual(int(result[0, -2]), 1)

    def cli_args(self, path, group):
        return SimpleNamespace(vcf=[path], region='chr1:1-11', upstream=0, downstream=0,
            type='fst', out=str(self.root / 'result'), het=None, maf=0, max_missing=1,
            group=str(group), size=5, step=5)

    def test_group_file_validation_and_literal_names(self):
        names = ('001', 'NA', 's3', 's4')
        _, path = self.variant([(1, 'A', 'G', DIFFERENT)], samples=names)
        group = self.root / 'groups.csv'
        args = self.cli_args(path, group)
        group.write_text('samples,groups\n001,01\nNA,01\ns3,NA\ns4,NA\n')
        view.run(args)
        frame = pd.read_csv(self.root / 'result.fst.csv', keep_default_na=False)
        self.assertEqual(frame.groups.iloc[0], '01_NA')
        self.close(float(frame.fst.iloc[0]), 1)
        invalid = ('samples,groups\n', 'samples,groups\n001,P\nNA,P\n',
                   'samples,groups\n001,P\n001,Q\n', 'samples,groups\n001,P\nmissing,Q\n',
                   'samples,groups\n001,P\n,Q\n', 'samples,groups\n001,P\ns3,\n',
                   'samples,groups,extra\n001,P,x\ns3,Q,x\n')
        for text in invalid:
            group.write_text(text)
            with self.subTest(text=text), self.assertRaises(ValueError):
                view.run(args)

    def test_unsupported_filters_fail_before_io(self):
        args = self.cli_args('nonexistent.vcf.gz', 'nonexistent.csv')
        for name, value in (('maf', .05), ('max_missing', .1)):
            setattr(args, name, value)
            with self.assertRaisesRegex(ValueError, 'not supported'):
                view.run(args)
            setattr(args, name, 0 if name == 'maf' else 1)

    def test_real_cli_dispatch_and_six_column_output(self):
        _, path = self.variant([(1, 'A', 'G', DIFFERENT), (11, 'A', 'G', HETERO)])
        group = self.root / 'groups.csv'
        group.write_text('samples,groups\ns1,P\ns2,P\ns3,Q\ns4,Q\n')
        tree = ast.parse(Path(find_spec('hastat.hastat').origin).read_text(encoding='utf-8'))
        main = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
        fake = SimpleNamespace(run=lambda args: None)
        namespace = dict(argparse=argparse, logging=logging, logger=logging.getLogger('fst_cli'),
                         view=view, stat=fake, plot=fake, gwas=fake, network=fake)
        exec(compile(ast.Module(body=[main], type_ignores=[]), 'hastat/hastat.py', 'exec'), namespace)
        before = set(self.root.iterdir())
        argv = ['hastat', 'view', '-v', path, '-r', 'chr1:1-11', '-g', str(group), '-t', 'fst',
                '--size', '5', '--step', '5', '-o', str(self.root / 'result')]
        with patch.object(sys, 'argv', argv):
            namespace['main']()
        self.assertEqual(set(self.root.iterdir()) - before, {self.root / 'result.fst.csv'})
        frame = pd.read_csv(self.root / 'result.fst.csv')
        self.assertEqual(list(frame.columns), ['chrom', 'start', 'end', 'counts', 'fst', 'groups'])
        self.close(frame.fst, [1, np.nan, 0])

    def test_multi_loci_api_and_cli(self):
        gv, path = self.variant([(1, 'A', 'G', DIFFERENT), (11, 'A', 'G', HETERO)])
        loci = [('chr1', 1, 5), ('chr1', 11, 11)]
        result = population.calculate_multi_loci_fst(gv, loci, ['s1', 's2'], ['s3', 's4'], size=5, step=5)
        expected = np.vstack([self.calculate(gv, left, right, 5, 5) for _, left, right in loci])
        np.testing.assert_array_equal(result, expected)
        group = self.root / 'groups.csv'
        group.write_text('samples,groups\ns1,P\ns2,P\ns3,Q\ns4,Q\n')
        annotation = self.root / 'genes.gff3'
        annotation.write_text('chr1\ttest\tgene\t1\t5\t.\t+\t.\tID=g1\nchr1\ttest\tgene\t11\t11\t.\t+\t.\tID=g2\n')
        args = self.cli_args(path, group)
        args.region, args.homo, args.gene_id, args.gff = None, 'g1,g2', None, str(annotation)
        view.run(args)
        frame = pd.read_csv(self.root / 'result.fst.csv')
        self.assertEqual(list(frame.columns), ['chrom', 'start', 'end', 'counts', 'fst', 'groups'])
        self.close(frame.fst, [1, 0])


if __name__ == '__main__':
    unittest.main(verbosity=2)
