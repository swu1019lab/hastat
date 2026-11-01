"""Small indexed VCF integration tests; execute inside a PBS allocation."""
import itertools
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import allel
import gffutils
import numpy as np
import pandas as pd
import pysam

from hastat.func import view
from hastat.utils.gene import GeneFeature, GeneVariant, TAJIMA_COLUMNS


COMPLETE = [
    (10, 'A', 'G', ['0/1', '0/1', '0/1']),
    (11, 'A', 'G', ['0/0', '0/1', '1/1']),
    (12, 'A', 'G', ['1/1', '1/1', '1/1']),
    (15, 'A', 'G', ['0/0', '0/0', '0/1']),
    (20, 'A', 'G', ['0/0', '0/1', '0/1']),
    (21, 'A', 'G', ['0/0', '0/0', '0/0']),
]


def write_vcf(path, rows, samples=('s1', 's2', 's3'), fmt='GT', chrom='chr1'):
    header = (f'##fileformat=VCFv4.2\n##contig=<ID={chrom},length=1000>\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples) + '\n')
    text = header + ''.join(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\t{fmt}\t' + '\t'.join(gt) + '\n'
                            for pos, ref, alt, gt in rows)
    path.write_text(text)
    pysam.tabix_compress(str(path), str(path) + '.gz', force=True)
    pysam.tabix_index(str(path) + '.gz', preset='vcf', force=True)
    return str(path) + '.gz'


def pairwise_pi(rows):
    """Independent enumeration of unordered chromosome pairs at each site."""
    total = 0.0
    for _, _, _, gt in rows:
        alleles = [int(a) for call in gt for a in call.replace('|', '/').split('/')]
        pairs = list(itertools.combinations(alleles, 2))
        total += sum(a != b for a, b in pairs) / len(pairs)
    return total


def tajima_reference(rows):
    """Tajima (1989) standardized theta difference; pair-enumerated numerator."""
    n = 2 * len(rows[0][3])
    s = sum(len(set('/'.join(gt).replace('|', '/').split('/'))) > 1 for _, _, _, gt in rows)
    if s < 3 or n < 4:
        return math.nan
    a1 = sum(1 / i for i in range(1, n))
    a2 = sum(1 / (i * i) for i in range(1, n))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    e1 = (b1 - 1 / a1) / a1
    e2 = (b2 - (n + 2) / (n * a1) + a2 / (a1 * a1)) / (a1 * a1 + a2)
    return (pairwise_pi(rows) - s / a1) / math.sqrt(e1 * s + e2 * s * (s - 1))


class PopulationValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.counter = 0

    def variant(self, rows=COMPLETE, **kwargs):
        self.counter += 1
        path = write_vcf(self.root / f'input{self.counter}.vcf', rows, **kwargs)
        handle = pysam.VariantFile(path)
        self.addCleanup(handle.close)
        return GeneVariant(handle)

    def test_both_boundaries_and_single_base(self):
        gv = self.variant()
        self.assertEqual(gv.get_geno_data('chr1', 10, 20).index.get_level_values('pos').tolist(), [10, 11, 12, 15, 20])
        self.assertEqual(gv.get_geno_data('chr1', 10, 10).index.get_level_values('pos').tolist(), [10])

    def test_fetch_overlapping_deletion_excluded_by_pos(self):
        gv = self.variant([(9, 'AAA', 'A', ['0/0', '0/1', '1/1']), COMPLETE[0]])
        self.assertEqual(gv.get_geno_data('chr1', 10, 10).index.get_level_values('pos').tolist(), [10])

    def test_gff_strand_extension_and_clipping(self):
        text = 'chr1\ttest\tgene\t10\t20\t.\t+\t.\tID=plus\nchr1\ttest\tgene\t10\t20\t.\t-\t.\tID=minus\n'
        db = gffutils.create_db(text, ':memory:', from_string=True)
        self.addCleanup(db.conn.close)
        gf = GeneFeature(db)
        self.assertEqual(gf.get_locus('plus', 3, 5), ('chr1', 7, 25))
        self.assertEqual(gf.get_locus('minus', 3, 5), ('chr1', 5, 23))
        self.assertEqual(gf.get_locus('plus', 50, 0), ('chr1', 1, 20))
        self.assertEqual(gf.get_locus('minus', 0, 50), ('chr1', 1, 20))
        with self.assertRaises(ValueError):
            gf.get_locus('plus', -1, 0)

    def test_partial_missing_symmetric_and_filters(self):
        calls = ['0/.', './0', '1/.', './1', './.', '.', '0/1', '1/1']
        gv = self.variant([(10, 'A', 'G', calls)], samples=tuple(f's{i}' for i in range(8)))
        self.assertEqual(gv.get_geno_data('chr1', 10, 10).iloc[0].tolist(), [-1] * 6 + [1, 2])
        self.assertEqual(gv.get_geno_data('chr1', 10, 10, 'code2').iloc[0].tolist(), [[-1, -1]] * 6 + [[0, 1], [1, 1]])
        for code in ('code1', 'code2'):
            self.assertEqual(len(gv.get_geno_data('chr1', 10, 10, code, maf_threshold=.25, max_missing_threshold=.75)), 1)
            self.assertTrue(gv.get_geno_data('chr1', 10, 10, code, maf_threshold=.251).empty)
            self.assertTrue(gv.get_geno_data('chr1', 10, 10, code, max_missing_threshold=.749).empty)
            self.assertTrue(gv.get_geno_data('chr1', 10, 10, code, het_threshold=.49).empty)

    def test_called_haploid_polyploid_and_missing_gt_rejected(self):
        for gt in ('0', '0/0/1', '././.'):
            gv = self.variant([(10, 'A', 'G', [gt, '0/0', '0/1'])])
            with self.assertRaisesRegex(ValueError, 'diploid GT'):
                gv.get_geno_data('chr1', 10, 10)
        gv = self.variant([(10, 'A', 'G', ['10', '10', '10'])], fmt='DP')
        with self.assertRaisesRegex(ValueError, 'Missing GT'):
            gv.get_geno_data('chr1', 10, 10)

    def test_empty_region_schema_and_windows(self):
        gv = self.variant()
        data = gv.get_geno_data('chr1', 50, 60, 'code2')
        self.assertEqual(data.shape, (0, 3))
        self.assertEqual(data.index.names, ['chrom', 'pos', 'ref', 'alt'])
        multi = gv.get_multi_loci_geno_data([('chr1', 50, 60)])
        self.assertEqual(multi.shape, (0, 3))
        pi = gv.get_pi_data('chr1', 50, 60, size=11, step=11)
        d = gv.get_tajima_d_data('chr1', 50, 60, size=11, step=11)
        np.testing.assert_allclose(pi[:, 1:].astype(float), [[50, 60, 11, 0, 0]])
        np.testing.assert_allclose(d[:, 1:].astype(float), [[50, 60, 0, np.nan]], equal_nan=True)

    def test_pi_all_diploid_genotype_combinations(self):
        for gt in itertools.product(('0/0', '0/1', '1/1'), repeat=3):
            rows = [(10, 'A', 'G', list(gt))]
            observed = float(self.variant(rows).get_pi_data('chr1', 10, 10)[0, -1])
            self.assertAlmostEqual(observed, pairwise_pi(rows), places=14)

    def test_pi_window_denominator_and_subset(self):
        gv = self.variant()
        pi = gv.get_pi_data('chr1', 10, 20, size=11, step=11)
        self.assertAlmostEqual(float(pi[0, -1]), pairwise_pi(COMPLETE[:5]) / 11, places=14)
        subset = [(p, r, a, [g[2], g[0]]) for p, r, a, g in COMPLETE[:5]]
        observed = gv.get_pi_data('chr1', 10, 20, size=11, step=11, sample_list=['s3', 's1'])
        self.assertAlmostEqual(float(observed[0, -1]), pairwise_pi(subset) / 11, places=14)

    def test_window_sizes_and_last_base(self):
        gv = self.variant()
        for method in (gv.get_pi_data, gv.get_tajima_d_data):
            result = method('chr1', 10, 20, size=5, step=5)
            np.testing.assert_array_equal(result[:, 1:3].astype(int), [[10, 14], [15, 19], [20, 20]])
            one_bp = method('chr1', 10, 12, size=1, step=1)
            np.testing.assert_array_equal(one_bp[:, 1:3].astype(int), [[10, 10], [11, 11], [12, 12]])
            overlap = method('chr1', 10, 20, size=5, step=3)
            np.testing.assert_array_equal(overlap[:, 1:3].astype(int), [[10, 14], [13, 17], [16, 20]])

    def test_phase_does_not_change_pi(self):
        gv = self.variant([(10, 'A', 'G', ['0|1', '1|0', '0/1'])])
        self.assertAlmostEqual(float(gv.get_pi_data('chr1', 10, 10)[0, -1]), .6)

    def test_missing_pi_uses_called_chromosomes_or_nan(self):
        gv = self.variant([(10, 'A', 'G', ['0/.', '0/0', '1/1']),
                           (20, 'A', 'G', ['./.', './0', '1/.'])])
        self.assertAlmostEqual(float(gv.get_pi_data('chr1', 10, 10)[0, -1]), 2 / 3)
        self.assertTrue(math.isnan(float(gv.get_pi_data('chr1', 20, 20)[0, -1])))

    def test_tajima_independent_formula_and_library_reference(self):
        rng = np.random.default_rng(20261031)
        for n in (2, 3, 5, 12):
            for _ in range(5):
                gt = rng.integers(0, 2, size=(12, n, 2))
                rows = [(i + 10, 'A', 'G', [f'{a}/{b}' for a, b in calls]) for i, calls in enumerate(gt)]
                gv = self.variant(rows, samples=tuple(f's{i}' for i in range(n)))
                observed = float(gv.get_tajima_d_data('chr1', 10, 30, size=21, step=21)[0, -1])
                np.testing.assert_allclose(observed, tajima_reference(rows), rtol=1e-12, atol=1e-12, equal_nan=True)
                reference = allel.windowed_tajima_d(np.arange(10, 22), allel.GenotypeArray(gt).count_alleles(), size=21, step=21, start=10, stop=30)[0][0]
                np.testing.assert_allclose(observed, reference, rtol=1e-12, atol=1e-12, equal_nan=True)

    def test_tajima_missing_filter_is_population_specific(self):
        rows = [(p, 'A', 'G', ['0/0', '0/1', '1/1']) for p in (10, 11, 12, 20, 21, 22)]
        rows[3] = (20, 'A', 'G', ['./0', '0/1', '1/1'])
        gv = self.variant(rows)
        data = gv.get_tajima_d_data('chr1', 10, 29, size=10, step=10)
        self.assertTrue(math.isfinite(float(data[0, -1])))
        self.assertTrue(math.isnan(float(data[1, -1])))
        subset = gv.get_tajima_d_data('chr1', 20, 29, size=10, step=10, sample_list=['s2', 's3'])
        self.assertTrue(math.isfinite(float(subset[0, -1])))
        expected = tajima_reference([(p, ref, alt, calls[1:]) for p, ref, alt, calls in rows[3:]])
        self.assertAlmostEqual(float(subset[0, -1]), expected, places=12)

    def test_tajima_excludes_incomplete_sites(self):
        complete = [(p, 'A', 'G', ['0/0', '0/1', '1/1']) for p in (10, 11, 12)]
        complete.append((13, 'A', 'G', ['1/1', '1/1', '1/1']))
        gv = self.variant(complete + [(14, 'A', 'G', ['0/.', '0/0', '1/1'])])
        with patch('hastat.utils.gene.allel.tajima_d', wraps=allel.tajima_d) as calculate:
            result = gv.get_tajima_d_data('chr1', 10, 19, size=10, step=10)
        self.assertEqual(result.shape, (1, 5))
        self.assertEqual(int(result[0, 3]), 5)
        self.assertAlmostEqual(float(result[0, -1]), tajima_reference(complete), places=12)
        np.testing.assert_array_equal(calculate.call_args.args[0], [[3, 3]] * 3 + [[0, 6]])

    def test_tajima_undefined_cases_return_nan(self):
        gv = self.variant([(10, 'A', 'G', ['./.', './.', './.']),
                           (20, 'A', 'G', ['0/0', '0/1', '1/1'])])
        for start, end, samples, count in [(50, 59, None, 0), (10, 19, None, 1),
                                           (20, 29, None, 1), (20, 29, ['s1'], 1)]:
            with self.subTest(start=start, samples=samples), patch('hastat.utils.gene.allel.tajima_d') as calculate:
                result = gv.get_tajima_d_data('chr1', start, end, size=10, step=10, sample_list=samples)
                self.assertEqual(result.shape, (1, 5))
                self.assertEqual(int(result[0, 3]), count)
                self.assertTrue(math.isnan(float(result[0, -1])))
                calculate.assert_not_called()

    def test_tajima_never_switches_sample_cohort_across_sites(self):
        rows = [(10, 'A', 'G', ['./.', '0/1', '1/1']),
                (11, 'A', 'G', ['0/0', './.', '1/1']),
                (12, 'A', 'G', ['0/0', '0/1', './.'])]
        gv = self.variant(rows)
        with patch('hastat.utils.gene.allel.tajima_d') as calculate:
            result = gv.get_tajima_d_data('chr1', 10, 19, size=10, step=10)
        self.assertEqual(int(result[0, 3]), 3)
        self.assertTrue(math.isnan(float(result[0, -1])))
        calculate.assert_not_called()

    def test_tajima_missing_matrices_match_independent_complete_site_reference(self):
        rng = np.random.default_rng(20261101)
        for n in (2, 5, 12):
            for _ in range(5):
                calls = rng.integers(0, 2, size=(20, n, 2))
                mask = rng.random((20, n)) < .08
                calls[mask] = -1
                rows = [(i + 10, 'A', 'G', ['./.' if a < 0 else f'{a}/{b}' for a, b in site])
                        for i, site in enumerate(calls)]
                keep = ~mask.any(axis=1)
                retained = [row for row, include in zip(rows, keep) if include]
                gv = self.variant(rows, samples=tuple(f's{i}' for i in range(n)))
                with patch('hastat.utils.gene.allel.tajima_d', wraps=allel.tajima_d) as calculate:
                    result = gv.get_tajima_d_data('chr1', 10, 29, size=20, step=20)
                expected = tajima_reference(retained) if retained else math.nan
                np.testing.assert_allclose(float(result[0, -1]), expected, rtol=1e-12, atol=1e-12, equal_nan=True)
                self.assertEqual(int(result[0, 3]), 20)
                if calculate.called:
                    expected_ac = np.column_stack([(calls[keep] == a).sum(axis=(1, 2)) for a in (0, 1)])
                    np.testing.assert_array_equal(calculate.call_args.args[0], expected_ac)

    def test_tajima_multiple_loci_and_overlapping_windows(self):
        rows = [(p, 'A', 'G', ['0/0', '0/1', '1/1']) for p in range(10, 20)]
        rows[4] = (14, 'A', 'G', ['0/.', '0/1', '1/1'])
        gv = self.variant(rows)
        result = gv.get_tajima_d_data('chr1', 10, 19, size=6, step=4)
        self.assertEqual(result.shape, (2, 5))
        np.testing.assert_array_equal(result[:, 3].astype(int), [6, 6])
        for row in result:
            retained = [r for r in rows if int(row[1]) <= r[0] <= int(row[2]) and r[0] != 14]
            self.assertAlmostEqual(float(row[-1]), tajima_reference(retained), places=12)
        loci = [('chr1', 10, 15), ('chr1', 14, 19)]
        np.testing.assert_array_equal(gv.get_multi_loci_tajima_d_data(loci, size=6, step=4), result)
        self.assertEqual(gv.get_multi_loci_tajima_d_data([]).shape, (0, len(TAJIMA_COLUMNS)))

    def test_tajima_small_sample_and_too_few_segregating_sites(self):
        gv = self.variant()
        self.assertTrue(math.isnan(float(gv.get_tajima_d_data('chr1', 10, 20, size=11, step=11, sample_list=['s1'])[0, -1])))
        self.assertTrue(math.isnan(float(gv.get_tajima_d_data('chr1', 10, 11, size=2, step=2)[0, -1])))

    def test_non_snp_excluded_and_duplicate_snp_rejected(self):
        rows = [(10, 'A', 'AT', ['0/0', '0/1', '1/1']), (11, 'A', 'G,T', ['0/0', '0/1', '1/2']), COMPLETE[2]]
        gv = self.variant(rows)
        self.assertEqual(len(gv.get_geno_data('chr1', 10, 12)), 2)
        result = gv.get_pi_data('chr1', 10, 12, size=3, step=3)
        self.assertEqual(float(result[0, -2]), 1)
        self.assertEqual(float(result[0, -1]), 0)
        duplicate = self.variant([COMPLETE[0], (10, 'A', 'T', ['0/0', '0/1', '1/1'])])
        with self.assertRaisesRegex(ValueError, 'one biallelic SNP'):
            duplicate.get_pi_data('chr1', 10, 10)

    def test_invalid_region_samples_and_windows(self):
        gv = self.variant()
        for start, end in ((0, 10), (11, 10), (1.5, 10)):
            with self.assertRaises(ValueError):
                gv.get_geno_data('chr1', start, end)
        for method in (gv.get_pi_data, gv.get_tajima_d_data):
            for sample_list in ([], ['absent'], ['s1', 's1']):
                with self.assertRaises(ValueError):
                    method('chr1', 10, 20, sample_list=sample_list)
            for size, step in ((0, 1), (1, 0), (1.5, 1)):
                with self.assertRaises(ValueError):
                    method('chr1', 10, 20, size=size, step=step)

    def test_multi_loci_matches_separate_regions(self):
        gv = self.variant()
        loci = [('chr1', 10, 12), ('chr1', 15, 21)]
        for single, multi in ((gv.get_pi_data, gv.get_multi_loci_pi_data), (gv.get_tajima_d_data, gv.get_multi_loci_tajima_d_data)):
            expected = np.vstack([single(*locus, size=10, step=10) for locus in loci])
            np.testing.assert_array_equal(multi(loci, size=10, step=10), expected)

    def test_fst_shared_reader_regression(self):
        rows = [(p, r, a, gt + ['0/1']) for p, r, a, gt in COMPLETE[:5]]
        gv = self.variant(rows, samples=('s1', 's2', 's3', 's4'))
        calls = np.array([[[int(a) for a in call.split('/')] for call in gt] for _, _, _, gt in rows])
        expected = allel.windowed_weir_cockerham_fst([r[0] for r in rows], allel.GenotypeArray(calls), [[0, 1], [2, 3]], size=11, step=11, start=10, stop=20)[0]
        original_err = np.geterr()
        try:
            actual = gv.get_fst_data('chr1', 10, 20, ['s1', 's2'], ['s3', 's4'], size=11, step=11)
        finally:
            np.seterr(**original_err)
        np.testing.assert_allclose(actual[:, -1].astype(float), expected, equal_nan=True)

    def test_view_region_and_population_csv_contract(self):
        path = write_vcf(self.root / 'cli.vcf', COMPLETE, chrom='chr-1')
        groups = self.root / 'groups.csv'
        groups.write_text('samples,groups\ns1,P\ns2,P\ns3,P\n')
        args = SimpleNamespace(vcf=[path], region='chr-1:11-19', upstream=1, downstream=1,
                               type='geno', out=str(self.root / 'result'), het=None, maf=0, max_missing=1,
                               group=str(groups), size=11, step=11)
        view.run(args)
        frame = pd.read_csv(self.root / 'result.geno.csv')
        self.assertEqual(frame.pos.tolist(), [10, 11, 12, 15, 20])
        for kind, cols in [('pi', ['chrom', 'start', 'end', 'n_bases', 'counts', 'pi', 'groups']),
                           ('tajima_d', TAJIMA_COLUMNS + ['groups'])]:
            args.type = kind
            view.run(args)
            frame = pd.read_csv(self.root / f'result.{kind}.csv')
            self.assertEqual(list(frame.columns), cols)
            expected = pairwise_pi(COMPLETE[:5]) / 11 if kind == 'pi' else tajima_reference(COMPLETE[:5])
            self.assertAlmostEqual(frame[kind].iloc[0], expected, places=12)
            if kind == 'tajima_d':
                self.assertNotIn('status', frame.columns)
                self.assertEqual(frame.counts.iloc[0], 5)
        # Missing-call handling must not add fields to the normal CSV.
        missing_rows = [(p, 'A', 'G', ['0/0', '0/1', '1/1']) for p in (10, 11, 12)]
        missing_rows.append((14, 'A', 'G', ['0/.', '0/1', '1/1']))
        args.vcf = [write_vcf(self.root / 'missing_cli.vcf', missing_rows)]
        args.region, args.upstream, args.downstream = 'chr1:10-19', 0, 0
        args.size, args.step = 10, 10
        view.run(args)
        frame = pd.read_csv(self.root / 'result.tajima_d.csv')
        self.assertEqual(list(frame.columns), TAJIMA_COLUMNS + ['groups'])
        self.assertEqual(frame.counts.iloc[0], 4)
        self.assertAlmostEqual(frame.tajima_d.iloc[0], tajima_reference(missing_rows[:3]), places=12)
        self.assertFalse(list(self.root.glob('*.summary.json')))
        args.maf = .05
        with self.assertRaisesRegex(ValueError, 'not supported'):
            view.run(args)


if __name__ == '__main__':
    unittest.main()
