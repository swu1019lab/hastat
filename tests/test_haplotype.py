"""Single/multi-locus genotype profiles and counts; execute inside PBS."""
import argparse
import ast
from importlib.util import find_spec
from collections import Counter
import logging
from pathlib import Path
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
NAMES = ('001', 'NA', 's3', 's4', 's5', 's6')
ROWS = [(10, 'A', 'G', ['0/0', '0/0', '0/0', '1/1', '1/1', '0/1']),
        (20, 'C', 'T', ['0/0', '0/0', '0/0', '1/1', '1/1', '0/1']),
        (30, 'G', 'A', ['0/1', '0/1', '0/1', '0/0', '0/0', '1/1'])]


def reference(rows, names):
    """Enumerate complete string-GT vectors and count using an independent Counter."""
    retained = []
    for pos, ref, alt, calls in sorted(rows):
        decoded = [call.replace('|', '/').split('/') for call in calls]
        if all(len(pair) == 2 and '.' not in pair for pair in decoded):
            retained.append((('chr1', pos, ref, alt), [sum(map(int, pair)) for pair in decoded]))
    keys = [key for key, _ in retained]
    if not retained:
        return keys, {}, {}
    profiles = {name: tuple(values[i] for _, values in retained) for i, name in enumerate(names)}
    counted = Counter(profiles.values())
    ranked = sorted(counted, key=lambda p: (-counted[p], tuple(-x for x in p)))
    labels = {p: f'Hap{i+1}' for i, p in enumerate(ranked)}
    return keys, {s: labels[p] for s, p in profiles.items()}, {labels[p]: p for p in ranked}


class HaplotypeValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root, self.counter = Path(self.tmp.name), 0

    def variant(self, rows=ROWS, names=NAMES):
        self.counter += 1
        path = write_vcf(self.root / f'input{self.counter}.vcf', rows, samples=names)
        handle = pysam.VariantFile(path)
        self.addCleanup(handle.close)
        return GeneVariant(handle), path

    def assert_reference(self, gv, rows, names):
        keys, assignments, profiles = reference(rows, names)
        table = gv.hap_table('chr1', 1, 100)
        groups = gv.hap_groups('chr1', 1, 100)
        self.assertEqual(list(table.iloc[:, :4].itertuples(index=False, name=None)), keys)
        self.assertEqual(list(table.columns[4:]), list(profiles))
        self.assertEqual(list(groups.columns), ['samples', 'haplotypes'])
        self.assertEqual(dict(zip(groups.samples, groups.haplotypes)), assignments)
        for hap, states in profiles.items():
            self.assertEqual(tuple(table[hap]), states)
        freq = gv.calc_hap_freq('chr1', 1, 100)
        self.assertEqual(list(freq.columns), ['samples'])
        self.assertEqual(freq.samples.to_dict(), dict(Counter(assignments.values())))
        self.assertTrue(all(np.issubdtype(dtype, np.integer) for dtype in freq.dtypes))
        loci = [('chr1', 1, 49), ('chr1', 50, 100)]
        pd.testing.assert_frame_equal(table, gv.multi_loci_hap_table(loci))
        pd.testing.assert_frame_equal(groups, gv.multi_loci_hap_groups(loci))
        pd.testing.assert_frame_equal(freq, gv.calc_multi_loci_hap_freq(loci))

    def metadata(self):
        return pd.DataFrame({'samples': NAMES, 'groups': ['P', 'P', 'Q', 'P', 'Q', 'Q']})

    def test_complete_table_group_and_counts(self):
        gv, _ = self.variant()
        self.assert_reference(gv, ROWS, NAMES)
        self.assertEqual(gv.calc_hap_freq('chr1', 1, 100).samples.tolist(), [3, 2, 1])

    def test_multi_region_grouped_counts_are_integer_counts(self):
        gv, _ = self.variant()
        expected = pd.DataFrame({'P': [2, 1, 0], 'Q': [1, 1, 1]}, index=pd.Index(['Hap1', 'Hap2', 'Hap3'], name='haplotypes'))
        single = gv.calc_hap_freq('chr1', 1, 100, self.metadata())
        multi = gv.calc_multi_loci_hap_freq([('chr1', 1, 10), ('chr1', 20, 30)], self.metadata())
        pd.testing.assert_frame_equal(single, expected)
        pd.testing.assert_frame_equal(multi, expected)
        self.assertEqual(single.sum().to_dict(), {'P': 3, 'Q': 3})

    def test_valid_subset_preserves_global_hap_ids(self):
        gv, _ = self.variant()
        metadata = pd.DataFrame({'samples': ['s6', 's4'], 'groups': ['Q', 'P']})
        result = gv.calc_hap_freq('chr1', 1, 100, metadata)
        self.assertEqual(list(result.index), ['Hap2', 'Hap3'])
        self.assertEqual(result.to_dict(), {'P': {'Hap2': 1, 'Hap3': 0}, 'Q': {'Hap2': 0, 'Hap3': 1}})
        self.assertEqual(gv.hap_table('chr1', 1, 100).shape[1], 7)

    def test_missing_does_not_create_profiles_or_change_raw_geno(self):
        rows = [(10, 'A', 'G', ['0/0'] * 6), (20, 'C', 'T', ['0/0', '0/.', '0/0', '0/0', '0/0', '0/0'])]
        gv, _ = self.variant(rows)
        before = gv.get_geno_data('chr1', 1, 100)
        self.assert_reference(gv, rows, NAMES)
        pd.testing.assert_frame_equal(before, gv.get_geno_data('chr1', 1, 100))
        self.assertEqual(before.loc[('chr1', 20, 'C', 'T'), 'NA'], -1)
        self.assertEqual(gv.calc_hap_freq('chr1', 1, 100).samples.to_dict(), {'Hap1': 6})
        # Even an unlisted sample's missing call excludes the locus: IDs match table/group.
        subset = self.metadata().iloc[[0, 2]]
        self.assertEqual(gv.calc_hap_freq('chr1', 1, 100, subset).sum(axis=1).to_dict(), {'Hap1': 2})

    def test_empty_and_no_complete_loci_have_no_phantom_hap(self):
        for rows in ([], [(10, 'A', 'G', ['./.'] * 6)], [(10, 'A', 'G', ['0/.'] + ['0/0'] * 5)]):
            gv, _ = self.variant(rows)
            self.assert_reference(gv, rows, NAMES)
            grouped = gv.calc_hap_freq('chr1', 1, 100, self.metadata())
            self.assertTrue(grouped.empty)
            self.assertEqual(list(grouped.columns), ['P', 'Q'])

    def test_heterozygosity_and_maf_filters_forwarded(self):
        rows = [(10, 'A', 'G', ['0/0'] * 5 + ['0/1']),
                (20, 'C', 'T', ['0/0'] * 3 + ['1/1'] * 3),
                (30, 'G', 'A', ['0/0'] * 5 + ['./.'])]
        gv, _ = self.variant(rows)
        for kwargs in ({'het_threshold': 0}, {'maf_threshold': .2}):
            table = gv.hap_table('chr1', 1, 100, **kwargs)
            self.assertEqual(table.pos.tolist(), [20])
            self.assertEqual(gv.calc_hap_freq('chr1', 1, 100, **kwargs).samples.tolist(), [3, 3])
            pd.testing.assert_frame_equal(table, gv.multi_loci_hap_table([('chr1', 1, 15), ('chr1', 20, 100)], **kwargs))
        self.assertEqual(gv.get_geno_data('chr1', 1, 100, max_missing_threshold=0).shape[0], 2)
        self.assertEqual(gv.hap_table('chr1', 1, 100, max_missing_threshold=0).pos.tolist(), [10, 20])

    def test_overlap_and_duplicate_variant_handling(self):
        gv, _ = self.variant()
        pd.testing.assert_frame_equal(gv.hap_table('chr1', 1, 100),
            gv.multi_loci_hap_table([('chr1', 20, 100), ('chr1', 1, 25), ('chr1', 10, 30)]))
        duplicate, _ = self.variant(ROWS[:1] + ROWS)
        with self.assertRaisesRegex(ValueError, 'Duplicate variant'):
            duplicate.hap_groups('chr1', 1, 100)
        with self.assertRaises(ValueError):
            gv.multi_loci_hap_table([])

    def test_order_and_ties_remain_consistent(self):
        calls = ['0/0', '0/1', '1/1'] * 2
        rows = [(10, 'A', 'G', calls), (30, 'C', 'T', ['1/1'] * 6)]
        gv, _ = self.variant(rows)
        other_rows = [(p, r, a, gt[::-1]) for p, r, a, gt in rows]
        shuffled, _ = self.variant(other_rows, NAMES[::-1])
        self.assert_reference(gv, rows, NAMES)
        self.assert_reference(shuffled, other_rows, NAMES[::-1])
        pd.testing.assert_frame_equal(gv.hap_table('chr1', 1, 100), shuffled.hap_table('chr1', 1, 100))
        self.assertEqual(tuple(gv.hap_table('chr1', 1, 100).Hap1), (2, 2))

    def test_more_than_nine_profiles_use_numeric_rank_order(self):
        names = tuple(f'x{i}' for i in range(12))
        states = [(i // 9, (i // 3) % 3, i % 3) for i in range(12)]
        labels = ['0/0', '0/1', '1/1']
        rows = [(pos, 'A', 'G', [labels[state[j]] for state in states]) for j, pos in enumerate((10, 20, 30))]
        gv, _ = self.variant(rows, names)
        self.assert_reference(gv, rows, names)
        self.assertEqual(list(gv.calc_hap_freq('chr1', 1, 100).index), [f'Hap{i}' for i in range(1, 13)])

    def test_invalid_metadata_is_not_silently_dropped_or_duplicated(self):
        gv, _ = self.variant()
        bad = [pd.DataFrame({'samples': ['001', '001'], 'groups': ['P', 'Q']}),
               pd.DataFrame({'samples': ['absent'], 'groups': ['P']}),
               pd.DataFrame({'samples': [''], 'groups': ['P']}),
               pd.DataFrame({'samples': ['001'], 'groups': ['']}),
               pd.DataFrame({'samples': ['001'], 'groups': ['haplotypes']}),
               pd.DataFrame(columns=['samples', 'groups']),
               pd.DataFrame({'samples': ['001'], 'groups': ['P'], 'extra': [1]})]
        for metadata in bad:
            with self.subTest(metadata=metadata), self.assertRaises(ValueError):
                gv.calc_hap_freq('chr1', 1, 100, metadata)

    def test_seeded_independent_profile_and_count_oracle(self):
        rng = np.random.default_rng(20261031)
        states = np.array(['0/0', '0/1', '1/1', './.'])
        for _ in range(30):
            rows = [(p, 'A', 'G', rng.choice(states, len(NAMES), p=[.3, .3, .35, .05]).tolist()) for p in range(10, 90, 10)]
            gv, _ = self.variant(rows)
            self.assert_reference(gv, rows, NAMES)
            _, assignment, _ = reference(rows, NAMES)
            actual = gv.calc_hap_freq('chr1', 1, 100, self.metadata())
            expected = Counter((assignment[s], g) for s, g in self.metadata().itertuples(index=False, name=None) if s in assignment)
            for hap in actual.index:
                for group in actual.columns:
                    self.assertEqual(actual.loc[hap, group], expected[hap, group])
            self.assertEqual(int(actual.to_numpy().sum()), len(assignment))

    def test_compare_uses_same_profile_builder(self):
        gv, _ = self.variant()
        other, _ = self.variant()
        compared = MultiPopulationAnalyzer([gv, other], ['P', 'Q']).compare_populations([('chr1', 1, 100)])
        single = gv.hap_table('chr1', 1, 100)
        merged = compared['merged_hap_table'].rename(columns=lambda name: name.replace('MergedHap', 'Hap'))
        pd.testing.assert_frame_equal(single, merged)

    def cli_args(self, path, groups=None):
        return SimpleNamespace(vcf=[path], region='chr1:1-100', upstream=0, downstream=0,
            type='table', out=str(self.root / 'result'), het=None, maf=0, max_missing=1, group=groups)

    def test_cli_flat_header_literal_ids_and_no_extra_files(self):
        _, path = self.variant()
        group = self.root / 'groups.csv'
        self.metadata().to_csv(group, index=False)
        args = self.cli_args(path, str(group))
        before = set(self.root.iterdir())
        for kind in ('table', 'group', 'freq'):
            args.type = kind
            view.run(args)
        self.assertEqual({p.name for p in set(self.root.iterdir()) - before}, {'result.table.csv', 'result.group.csv', 'result.freq.csv'})
        freq = pd.read_csv(self.root / 'result.freq.csv')
        self.assertEqual(list(freq.columns), ['haplotypes', 'P', 'Q'])
        self.assertEqual(freq.P.tolist(), [2, 1, 0])
        self.assertEqual(freq.Q.tolist(), [1, 1, 1])
        assignments = pd.read_csv(self.root / 'result.group.csv', keep_default_na=False, dtype=str)
        self.assertEqual(assignments.samples.tolist(), list(NAMES))
        network.run(SimpleNamespace(hap_table=str(self.root / 'result.table.csv'), method='MSN', out=str(self.root / 'network.csv')))
        self.assertEqual(read_network(self.root / 'network.csv').number_of_nodes(), 3)

    def test_multi_loci_cli_actual_parser_and_grouped_freq(self):
        _, path = self.variant()
        group = self.root / 'groups.csv'
        self.metadata().to_csv(group, index=False)
        gff = self.root / 'genes.gff3'
        gff.write_text('chr1\ttest\tgene\t10\t10\t.\t+\t.\tID=g1\nchr1\ttest\tgene\t30\t30\t.\t+\t.\tID=g2\n')
        tree = ast.parse(Path(find_spec('hastat.hastat').origin).read_text(encoding='utf-8'))
        main = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
        fake = SimpleNamespace(run=lambda args: None)
        namespace = dict(argparse=argparse, logging=logging, logger=logging.getLogger('hap_cli'),
                         view=view, stat=fake, plot=fake, gwas=fake, network=fake)
        exec(compile(ast.Module(body=[main], type_ignores=[]), 'hastat/hastat.py', 'exec'), namespace)
        argv = ['hastat', 'view', '-v', path, '-a', str(gff), '--homo', 'g1,g2', '-t', 'freq',
                '-g', str(group), '-o', str(self.root / 'multi')]
        with patch.object(sys, 'argv', argv):
            namespace['main']()
        frame = pd.read_csv(self.root / 'multi.freq.csv')
        self.assertEqual(list(frame.columns), ['haplotypes', 'P', 'Q'])
        self.assertEqual(frame[['P', 'Q']].sum().to_dict(), {'P': 3, 'Q': 3})

    def test_empty_reruns_replace_stale_outputs(self):
        _, path = self.variant()
        _, empty = self.variant([])
        group = self.root / 'groups.csv'
        self.metadata().to_csv(group, index=False)
        args = self.cli_args(path, str(group))
        for source in (path, empty):
            args.vcf = [source]
            for kind in ('table', 'group', 'freq'):
                args.type = kind
                view.run(args)
        expected = {'table': ['chrom', 'pos', 'ref', 'alt'], 'group': ['samples', 'haplotypes'], 'freq': ['haplotypes', 'P', 'Q']}
        for kind, columns in expected.items():
            frame = pd.read_csv(self.root / f'result.{kind}.csv')
            self.assertTrue(frame.empty)
            self.assertEqual(list(frame.columns), columns)
        self.assertFalse(list(self.root.glob('*.summary.json')))

    def test_multiple_chromosomes_and_reordered_loci(self):
        path = self.root / 'two.vcf'
        path.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##contig=<ID=chr2,length=100>\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n'
            'chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\nchr2\t20\t.\tC\tT\t.\tPASS\t.\tGT\t1/1\t0/0\n')
        pysam.tabix_compress(str(path), str(path) + '.gz', force=True)
        pysam.tabix_index(str(path) + '.gz', preset='vcf', force=True)
        handle = pysam.VariantFile(str(path) + '.gz')
        self.addCleanup(handle.close)
        gv = GeneVariant(handle)
        loci = [('chr2', 1, 50), ('chr1', 1, 50)]
        table = gv.multi_loci_hap_table(loci)
        self.assertEqual(table.chrom.tolist(), ['chr1', 'chr2'])
        self.assertEqual(tuple(table.Hap1), (2, 0))
        self.assertEqual(gv.calc_multi_loci_hap_freq(loci).samples.to_dict(), {'Hap1': 1, 'Hap2': 1})


if __name__ == '__main__':
    unittest.main(verbosity=2)
