"""Independent profile oracle and real cross-VCF integration; run in PBS."""
import argparse
import ast
from importlib.util import find_spec
from collections import Counter
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
import pysam

from hastat.func import view, network
from hastat.utils.gene import GeneVariant, MultiPopulationAnalyzer
from hastat.utils.network_io import read_network
from test_population import write_vcf

ROOT = Path(__file__).resolve().parents[1]
META = ['chrom', 'pos', 'ref', 'alt']
NAMES = ('s1', 's2', 's3')
P = [(10, 'A', 'G', ['0/0', '0/0', '1/1']), (20, 'C', 'T', ['0/1', '0/1', '1/1']),
     (30, 'G', 'A', ['1/1', '1/1', '0/0'])]
Q = [(10, 'A', 'G', ['1/1', '0/0', '0/0']), (20, 'C', 'T', ['1/1', '0/1', '0/1']),
     (30, 'G', 'A', ['0/0', '1/1', '1/1'])]


def oracle(populations):
    """Use dictionaries and complete string GT calls; no production grouping."""
    record_maps = []
    for rows, samples in populations:
        records = {}
        for pos, ref, alt, calls in rows:
            key = ('chr1', pos, ref.upper() if set(ref.upper()) <= set('ACGTN') else ref,
                   alt.upper() if set(alt.upper()) <= set('ACGTN') else alt)
            values = []
            for call in calls:
                alleles = call.replace('|', '/').split('/')
                values.append(sum(map(int, alleles)) if len(alleles) == 2 and '.' not in alleles else None)
            records[key] = values
        record_maps.append(records)
    shared = set.intersection(*(set(records) for records in record_maps))
    usable = sorted(key for key in shared if all(None not in records[key] for records in record_maps))
    profiles = []
    for i, ((_, samples), records) in enumerate(zip(populations, record_maps), 1):
        profiles.append({f'{sample}_{i}': tuple(records[key][j] for key in usable)
                         for j, sample in enumerate(samples)} if usable else {})
    return usable, profiles


def check_result(test, result, populations, names):
    keys, profiles = oracle(populations)
    merged = {sample: profile for part in profiles for sample, profile in part.items()}
    for label, prefix, values in [(name, name + 'Hap', part) for name, part in zip(names, profiles)] + [('merged', 'MergedHap', merged)]:
        table, groups = result[f'{label}_hap_table'], result[f'{label}_hap_groups']
        test.assertEqual(list(table.columns[:4]), META)
        test.assertEqual(list(groups.columns), ['samples', 'haplotypes'])
        test.assertEqual(list(table.loc[:, META].itertuples(index=False, name=None)), keys)
        counts = Counter(values.values())
        ranked = sorted(counts, key=lambda p: (-counts[p], tuple(-x for x in p)))
        labels = {profile: f'{prefix}{i+1}' for i, profile in enumerate(ranked)}
        test.assertEqual(dict(zip(groups.samples, groups.haplotypes)), {s: labels[p] for s, p in values.items()})
        test.assertEqual(list(table.columns[4:]), [f'{prefix}{i+1}' for i in range(len(ranked))])
        for profile, hap in labels.items():
            test.assertEqual(tuple(table[hap]), profile)


class CompareValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root, self.counter = Path(self.tmp.name), 0

    def variants(self, populations):
        variants, paths = [], []
        for rows, samples in populations:
            self.counter += 1
            path = write_vcf(self.root / f'input{self.counter}.vcf', rows, samples=samples)
            handle = pysam.VariantFile(path)
            self.addCleanup(handle.close)
            variants.append(GeneVariant(handle))
            paths.append(path)
        return variants, paths

    def analyze(self, populations, names=('P', 'Q'), loci=(('chr1', 1, 100),)):
        variants, _ = self.variants(populations)
        return MultiPopulationAnalyzer(variants, names).compare_populations(loci)

    def test_complete_profiles_and_merged_identity(self):
        populations = [(P, NAMES), (Q, NAMES)]
        result = self.analyze(populations)
        check_result(self, result, populations, ('P', 'Q'))
        self.assertEqual(result['merged_hap_table'].shape[1] - 4, 2)
        self.assertEqual(result['merged_hap_groups'].haplotypes.value_counts().tolist(), [4, 2])

    def test_missing_calls_do_not_form_extra_haplotypes(self):
        populations = [([(10, 'A', 'G', ['0/0'] * 3), (20, 'C', 'T', ['0/0', './.', '0/0'])], NAMES),
                       ([(10, 'A', 'G', ['0/0'] * 3), (20, 'C', 'T', ['0/0'] * 3)], NAMES)]
        result = self.analyze(populations)
        check_result(self, result, populations, ('P', 'Q'))
        self.assertEqual(result['merged_hap_table'].pos.tolist(), [10])
        self.assertEqual(result['merged_hap_groups'].haplotypes.nunique(), 1)
        self.assertEqual(result['variant_class'].present_in_all.sum(), 2)

    def test_incomplete_shared_loci_produce_no_assignments(self):
        populations = [([(10, 'A', 'G', ['0/.', '0/0', '1/1'])], NAMES),
                       ([(10, 'A', 'G', ['0/0'] * 3)], NAMES)]
        result = self.analyze(populations)
        check_result(self, result, populations, ('P', 'Q'))
        self.assertTrue(result['merged_hap_groups'].empty)

    def test_one_empty_population_retains_samples_and_presence_column(self):
        populations = [(P, NAMES), ([], ('e1', 'e2'))]
        result = self.analyze(populations)
        check_result(self, result, populations, ('P', 'Q'))
        self.assertFalse(result['variant_class'].present_in_Q.any())
        metrics = result['comparison_summary'].set_index('metric').value
        self.assertEqual(metrics['Q samples'], 2)
        self.assertEqual(metrics['Q total variants'], 0)

    def test_all_empty_and_no_shared_records(self):
        for populations in ([( [], NAMES), ([], NAMES)], [(P[:1], NAMES), (Q[1:], NAMES)]):
            result = self.analyze(populations)
            check_result(self, result, populations, ('P', 'Q'))
            self.assertEqual(list(result['variant_class'].columns), META + ['present_in_P', 'present_in_Q',
                'present_in_populations', 'num_populations', 'present_in_all', 'population_specific', 'specific_to'])
            self.assertEqual(len(result['comparison_summary']), 9)

    def test_exact_alt_matching_and_no_absent_reference_imputation(self):
        populations = [(P[:2], NAMES), ([Q[0], (20, 'C', 'G', ['0/0'] * 3)], NAMES)]
        result = self.analyze(populations)
        check_result(self, result, populations, ('P', 'Q'))
        self.assertEqual(result['merged_hap_table'].pos.tolist(), [10])
        self.assertEqual(result['variant_class'].population_specific.sum(), 2)

    def test_sequence_case_normalized_symbolic_ids_preserved(self):
        p = [(10, 'a', 'g', ['0/0', '0/1', '1/1']), (20, 'A', '<DEL>', ['0/0'] * 3)]
        q = [(10, 'A', 'G', ['0/0', '0/1', '1/1']), (20, 'A', '<del>', ['0/0'] * 3)]
        populations = [(p, NAMES), (q, NAMES)]
        result = self.analyze(populations)
        check_result(self, result, populations, ('P', 'Q'))
        self.assertEqual(result['merged_hap_table'].pos.tolist(), [10])

    def test_reference_swap_and_duplicate_records_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Incompatible SNP REF'):
            self.analyze([(P[:1], NAMES), ([(10, 'G', 'A', ['1/1'] * 3)], NAMES)])
        with self.assertRaisesRegex(ValueError, 'Duplicate variant'):
            self.analyze([(P[:1] + P[:1], NAMES), (Q, NAMES)])

    def test_overlapping_loci_are_fetched_once(self):
        populations = [(P, NAMES), (Q, NAMES)]
        result = self.analyze(populations, loci=[('chr1', 20, 50), ('chr1', 1, 25), ('chr1', 1, 10)])
        check_result(self, result, populations, ('P', 'Q'))
        self.assertEqual(result['merged_hap_table'].pos.tolist(), [10, 20, 30])

    def test_profile_ids_stable_under_sample_and_vcf_order(self):
        populations = [(P, NAMES), (Q, NAMES)]
        forward = self.analyze(populations)
        shuffled = [([(p, r, a, [gt[2], gt[0], gt[1]]) for p, r, a, gt in rows], ('s3', 's1', 's2'))
                    for rows, _ in populations]
        result = self.analyze(shuffled)
        check_result(self, result, shuffled, ('P', 'Q'))
        pd.testing.assert_frame_equal(forward['merged_hap_table'], result['merged_hap_table'])
        reverse = self.analyze(populations[::-1], names=('Q', 'P'))
        pd.testing.assert_frame_equal(forward['merged_hap_table'], reverse['merged_hap_table'])
        check_result(self, reverse, populations[::-1], ('Q', 'P'))

    def test_tie_breaking_and_row_order(self):
        idx = pd.MultiIndex.from_tuples([('chr1', 30, 'A', 'G'), ('chr1', 10, 'A', 'G')], names=META)
        data = pd.DataFrame([[0, 1, 2, 0], [2, 0, 1, 0]], index=idx, columns=list('abcd'))
        class Stub:
            def __init__(self, frame):
                self.frame = frame
            def get_geno_data(self, *args):
                return self.frame
        first = MultiPopulationAnalyzer([Stub(data), Stub(data)], ['P', 'Q']).compare_populations([('chr1', 1, 40)])
        second = MultiPopulationAnalyzer([Stub(data.iloc[::-1, ::-1]), Stub(data.iloc[::-1])], ['P', 'Q']).compare_populations([('chr1', 1, 40)])
        pd.testing.assert_frame_equal(first['merged_hap_table'], second['merged_hap_table'])
        self.assertEqual(tuple(first['merged_hap_table'].MergedHap1), (2, 0))

    def test_seeded_independent_oracle(self):
        rng = np.random.default_rng(20261031)
        states = np.array(['0/0', '0/1', '1/1', './.'])
        for _ in range(30):
            populations = []
            count = int(rng.integers(2, 5))
            labels = [f'P{i+1}' for i in range(count)]
            for pop in range(count):
                names = tuple(f'x{i}' for i in range(int(rng.integers(2, 6))))
                rows = [(pos, 'A', 'G', rng.choice(states, len(names), p=[.35, .3, .3, .05]).tolist())
                        for pos in range(10, 70, 10) if pos == 10 or rng.random() > .1]
                populations.append((rows, names))
            check_result(self, self.analyze(populations, names=labels), populations, labels)

    def test_population_names_and_loci_validation(self):
        variants, _ = self.variants([(P, NAMES), (Q, NAMES)])
        for labels in (['P'], ['P', 'P'], ['', 'Q'], ['merged', 'Q'], ['../P', 'Q'], ['P ', 'Q']):
            with self.subTest(labels=labels), self.assertRaises(ValueError):
                MultiPopulationAnalyzer(variants, labels)
        analyzer = MultiPopulationAnalyzer(variants, ['群体P', 'Q'])
        for loci in ([], [('chr1', 0, 10)], [('chr1', 20, 10)]):
            with self.assertRaises(ValueError):
                analyzer.compare_populations(loci)

    def cli_args(self, paths):
        return SimpleNamespace(vcf=paths, pop_names=['P', 'Q'], region='chr1:1-100', upstream=0,
            downstream=0, type='compare', out=str(self.root / 'result'), het=None, maf=0, max_missing=1)

    def test_unsupported_filters_rejected_before_io(self):
        args = self.cli_args(['absent1', 'absent2'])
        for key, value, default in (('het', .1, None), ('maf', .05, 0), ('max_missing', .2, 1)):
            setattr(args, key, value)
            with self.assertRaisesRegex(ValueError, 'not supported'):
                view.run(args)
            setattr(args, key, default)

    def test_real_cli_and_output_contract(self):
        _, paths = self.variants([(P, NAMES), (Q, NAMES)])
        tree = ast.parse(Path(find_spec('hastat.hastat').origin).read_text(encoding='utf-8'))
        main = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
        fake = SimpleNamespace(run=lambda args: None)
        namespace = dict(argparse=argparse, logging=logging, logger=logging.getLogger('compare_cli'),
                         view=view, stat=fake, plot=fake, gwas=fake, network=fake)
        exec(compile(ast.Module(body=[main], type_ignores=[]), 'hastat/hastat.py', 'exec'), namespace)
        argv = ['hastat', 'view', '-t', 'compare', '-v'] + paths + ['-n', 'P', 'Q', '-r', 'chr1:1-100', '-o', str(self.root / 'result')]
        existing = set(self.root.iterdir())
        with patch.object(sys, 'argv', argv):
            namespace['main']()
        names = {p.name for p in set(self.root.iterdir()) - existing}
        self.assertEqual(names, {'result_P.table.csv', 'result_P.group.csv', 'result_Q.table.csv', 'result_Q.group.csv',
            'result.merged.table.csv', 'result.merged.group.csv', 'result.variant_class.csv', 'result.comparison_summary.csv'})
        group = pd.read_csv(self.root / 'result.merged.group.csv')
        self.assertEqual(list(group.columns), ['samples', 'haplotypes'])
        self.assertEqual(set(group.samples), {f'{s}_{i}' for i in (1, 2) for s in NAMES})
        network.run(SimpleNamespace(hap_table=str(self.root / 'result.merged.table.csv'), method='MSN',
                                    out=str(self.root / 'network.csv')))
        self.assertEqual(read_network(self.root / 'network.csv').number_of_nodes(), 2)

    def test_empty_rerun_replaces_stale_assignments(self):
        _, paths = self.variants([(P, NAMES), (Q, NAMES)])
        args = self.cli_args(paths)
        view.run(args)
        _, args.vcf = self.variants([([], NAMES), ([], NAMES)])
        view.run(args)
        for label in ('_P', '_Q', '.merged'):
            table = pd.read_csv(self.root / f'result{label}.table.csv')
            groups = pd.read_csv(self.root / f'result{label}.group.csv')
            self.assertTrue(table.empty and groups.empty)
            self.assertEqual(list(table.columns), META)
            self.assertEqual(list(groups.columns), ['samples', 'haplotypes'])
        self.assertEqual(len(list(self.root.glob('result*.csv'))), 8)
        self.assertFalse(list(self.root.glob('*.summary.json')))

    def test_hash_seed_independence(self):
        _, paths = self.variants([(P, NAMES), (Q, NAMES)])
        script = ('import json,pysam; from hastat.utils.gene import GeneVariant,MultiPopulationAnalyzer; '
                  f'v=[GeneVariant(pysam.VariantFile(p)) for p in {paths!r}]; '
                  'r=MultiPopulationAnalyzer(v,["P","Q"]).compare_populations([("chr1",1,100)]); '
                  'print(json.dumps({k:x.to_csv(index=False) for k,x in r.items()},sort_keys=True))')
        runs = [subprocess.run([sys.executable, '-c', script], cwd=ROOT,
                env=dict(os.environ, PYTHONHASHSEED=seed), check=True, capture_output=True, text=True).stdout
                for seed in ('1', '7', '20261031')]
        self.assertEqual(len(set(runs)), 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
