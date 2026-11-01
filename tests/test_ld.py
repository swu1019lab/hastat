"""LD boundary/VCF/render checks and independent centered-dot-product oracle; PBS only."""
from itertools import product
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import allel
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

from hastat.utils.gene import GeneVariant
from hastat.utils.haplotype import normalize_loci
from hastat.utils.ld import ld_r_squared
from hastat.viz import ldheatmap
from test_population import write_vcf

DOSAGE = np.array([[0, 1, 2, 0], [0, 1, 2, 0], [2, 1, 0, 2], [0, 0, 1, 2]], dtype=np.int8)


def reference(data):
    """Pearson r² using centered vectors, independent of the library's moment sums."""
    data = np.asarray(data, dtype=float)
    matrix = np.full((len(data), len(data)), np.nan)
    for i in range(len(data)):
        for j in range(i, len(data)):
            keep = np.isfinite(data[i]) & np.isfinite(data[j]) & (data[i] >= 0) & (data[j] >= 0)
            if keep.sum() < 2:
                continue
            x, y = data[i, keep], data[j, keep]
            x, y = x - x.mean(), y - y.mean()
            denominator = np.dot(x, x) * np.dot(y, y)
            if denominator > 0:
                matrix[i, j] = matrix[j, i] = np.dot(x, y)**2 / denominator
    return matrix


def rows_from_dosage(data, first=10):
    calls = {-1: './.', 0: '0/0', 1: '0/1', 2: '1/1'}
    return [(first + i*10, 'A', 'G', [calls[int(x)] for x in row]) for i, row in enumerate(data)]


class LDValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.addCleanup(plt.close, 'all')
        self.root = Path(self.tmp.name)
        self.serial = 0

    def variant(self, rows=None, samples=('s1', 's2', 's3', 's4')):
        self.serial += 1
        path = write_vcf(self.root / f'input{self.serial}.vcf', rows_from_dosage(DOSAGE) if rows is None else rows,
                         samples=samples)
        handle = pysam.VariantFile(path)
        self.addCleanup(handle.close)
        return GeneVariant(handle)

    def close(self, actual, expected):
        np.testing.assert_allclose(actual, expected, rtol=2e-5, atol=2e-6, equal_nan=True)

    def test_complete_correlated_inverse_and_nonperfect_pairs(self):
        result = ld_r_squared(DOSAGE)
        self.close(result, reference(DOSAGE))
        self.assertEqual(result.shape, (4, 4))
        np.testing.assert_array_equal(np.diag(result), 1)
        self.assertEqual(result[0, 1], 1)
        self.assertEqual(result[0, 2], 1)
        self.assertLess(result[0, 3], 1)

    def test_all_four_sample_dosage_and_missing_patterns(self):
        data = np.array(list(product((-1, 0, 1, 2), repeat=4)), dtype=np.int8)
        result = ld_r_squared(data)
        self.close(result, reference(data))
        self.assertEqual(result.shape, (256, 256))
        self.assertTrue(np.all((result[np.isfinite(result)] >= 0) & (result[np.isfinite(result)] <= 1)))

    def test_empty_matrix_has_zero_rows_and_columns(self):
        for shape in ((0, 0), (0, 4)):
            self.assertEqual(ld_r_squared(np.empty(shape)).shape, (0, 0))

    def test_zero_or_one_sample_has_no_estimable_entries(self):
        for data in (np.empty((3, 0)), [[0], [1], [2]]):
            result = ld_r_squared(data)
            self.assertEqual(result.shape, (3, 3))
            self.assertTrue(np.isnan(result).all())

    def test_single_variant_diagonal_depends_on_dosage_variation(self):
        for calls, expected in (([0, 1, -1], 1), ([1, 1, 1], np.nan), ([2, 2, 2], np.nan),
                                ([0, 0, 0], np.nan), ([-1, -1, -1], np.nan), ([1, -1, -1], np.nan)):
            result = ld_r_squared([calls])
            self.assertEqual(result.shape, (1, 1))
            self.close(result, [[expected]])

    def test_two_variants_scalar_library_output_normalized(self):
        self.assertEqual(np.asarray(allel.rogers_huff_r(DOSAGE[:2])).ndim, 0)
        self.close(ld_r_squared(DOSAGE[:2]), [[1, 1], [1, 1]])

    def test_pairwise_missing_not_imputed_or_globally_deleted(self):
        data = np.array([[0, 2, 0, 2, -1, -1], [0, 2, -1, -1, 0, 2], [-1, -1, 0, 2, 0, 2]])
        # No globally complete sample, but each pair shares two variable calls.
        self.assertFalse(np.any(np.all(data >= 0, axis=0)))
        self.close(ld_r_squared(data), reference(data))
        self.assertTrue(np.isfinite(ld_r_squared(data)).all())

    def test_pair_constant_on_overlap_is_nan_despite_valid_diagonals(self):
        data = [[0, 0, 2, -1], [1, 2, -1, 0]]
        self.close(ld_r_squared(data), [[1, np.nan], [np.nan, 1]])
        self.close(ld_r_squared([[0, 1, -1], [-1, 1, 2]]), [[1, np.nan], [np.nan, 1]])

    def test_nan_input_equivalent_to_missing_code(self):
        data = DOSAGE.astype(float)
        data[0, 1] = np.nan
        changed = data.copy()
        changed[np.isnan(changed)] = -1
        self.close(ld_r_squared(data), ld_r_squared(changed))

    def test_invalid_states_or_shape_rejected(self):
        for data in ([0, 1, 2], [[0, 3]], [[0, .5]], [[0, -2]], [[0, np.inf]], [['x', 'y']]):
            with self.assertRaises(ValueError):
                ld_r_squared(data)

    def test_real_library_used_and_nonfinite_not_clipped_to_one(self):
        original = allel.rogers_huff_r
        with patch('hastat.utils.ld.allel.rogers_huff_r', wraps=original) as call:
            ld_r_squared(DOSAGE)
            call.assert_called_once()
            self.assertEqual(call.call_args.args[0].dtype, np.int8)
        with patch('hastat.utils.ld.allel.rogers_huff_r', return_value=np.array([1.0000001, np.inf, np.nan], dtype=np.float32)):
            result = ld_r_squared(DOSAGE[:3])
            self.assertEqual(result[0, 1], 1)
            self.assertTrue(np.isnan(result[0, 2]))
            self.assertTrue(np.isnan(result[1, 2]))

    def test_single_region_boundaries_empty_and_partial_calls(self):
        gv = self.variant()
        self.close(gv.get_ld_data('chr1', 10, 40), reference(DOSAGE))
        self.close(gv.get_ld_data('chr1', 10, 10), [[1]])
        self.assertEqual(gv.get_ld_data('chr1', 41, 100).shape, (0, 0))
        rows = [(10, 'A', 'G', ['0/.', '0/0', '0/1', '1/1']),
                (20, 'A', 'G', ['./1', '0/0', '0/1', '1/1'])]
        self.close(self.variant(rows).get_ld_data('chr1', 1, 30), [[1, 1], [1, 1]])

    def test_overlapping_reordered_multi_regions_count_once(self):
        gv = self.variant()
        loci = [('chr1', 20, 40), ('chr1', 1, 30), ('chr1', 10, 10)]
        self.close(gv.get_multi_loci_ld_data(loci), gv.get_ld_data('chr1', 1, 40))
        self.close(gv.get_multi_loci_ld_data(reversed(loci)), reference(DOSAGE))
        self.assertEqual(gv.get_multi_loci_ld_data([('chr1', 60, 80), ('chr1', 90, 100)]).shape, (0, 0))
        self.close(gv.get_multi_loci_ld_data([('chr1', 10, 10), ('chr1', 90, 100)]), [[1]])

    def test_empty_or_invalid_region_requests_rejected(self):
        gv = self.variant()
        for loci in ([], [('chr1', 0, 10)], [('chr1', 20, 10)], [('chr1', 1.5, 10)]):
            with self.assertRaises(ValueError):
                gv.get_multi_loci_ld_data(loci)

    def test_duplicate_records_rejected_but_indels_retained(self):
        row = rows_from_dosage(DOSAGE[:1])[0]
        gv = self.variant([row, row])
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            gv.get_ld_data('chr1', 1, 50)
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            gv.get_multi_loci_ld_data([('chr1', 1, 50)])
        self.close(self.variant([(10, 'AT', 'A', ['0/0', '0/1', '1/1', './.'])]).get_ld_data('chr1', 10, 10), [[1]])

    def test_cross_chromosome_matrix_order(self):
        path = self.root / 'multi.vcf'
        rows = rows_from_dosage(DOSAGE[:2])
        write_vcf(path, rows, samples=('s1', 's2', 's3', 's4'))
        text = path.read_text().replace('##FORMAT=<ID=GT', '##contig=<ID=chr2,length=1000>\n##FORMAT=<ID=GT')
        text += 'chr2\t10\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0\t1/1\n'
        path.write_text(text)
        pysam.tabix_compress(str(path), str(path)+'.gz', force=True)
        pysam.tabix_index(str(path)+'.gz', preset='vcf', force=True)
        with pysam.VariantFile(str(path)+'.gz') as handle:
            gv = GeneVariant(handle)
            loci = [('chr2', 1, 40), ('chr1', 1, 40)]
            frame = gv.get_multi_loci_geno_data(normalize_loci(loci))
            self.assertEqual(list(frame.index.get_level_values('chrom')), ['chr1', 'chr1', 'chr2'])
            self.close(gv.get_multi_loci_ld_data(loci), reference(frame.to_numpy()))

    def test_heatmap_masks_nan_and_skips_numeric_annotation(self):
        matrix = ld_r_squared([[0, 1, 2], [1, 1, 1], [2, 1, 0]])
        fig, ax = plt.subplots()
        ldheatmap.plot(pd.DataFrame(matrix), plot_value=True, ax=ax)
        fig.canvas.draw()
        self.assertEqual(int(np.ma.getmaskarray(ax.collections[0].get_array()).sum()), 3)
        self.assertEqual(len(ax.texts), 3)
        self.assertEqual([t.get_text() for t in ax.get_legend().get_texts()], ['Not estimable'])
        self.assertTrue(all(t.get_text() == '1.00' for t in ax.texts))
        fig.savefig(self.root / 'masked.pdf')

    def test_heatmap_single_position_and_repeated_positions_are_finite(self):
        for positions, matrix in (([10], [[1]]), ([10], [[np.nan]]), ([10, 10], [[1, 1], [1, 1]])):
            if len(positions) == 1:
                index = pd.Index(positions, name='pos')
            else:
                index = pd.MultiIndex.from_tuples([('chr1', 10, 'A', 'G'), ('chr1', 10, 'A', 'T')],
                                                names=['chrom', 'pos', 'ref', 'alt'])
            fig, ax = plt.subplots()
            track = ldheatmap.plot(pd.DataFrame(matrix, index=index, columns=index), plot_value=True, plot_snp=True, ax=ax)
            fig.canvas.draw()
            self.assertTrue(np.isfinite(track.collections[0].get_segments()).all())

    def test_heatmap_all_nan_and_hidden_diagonal(self):
        fig, ax = plt.subplots()
        ldheatmap.plot(pd.DataFrame(np.full((2, 2), np.nan)), plot_value=True, plot_diag=False, ax=ax)
        fig.canvas.draw()
        self.assertEqual(len(ax.collections[0].get_array()), 1)
        self.assertEqual(len(ax.texts), 0)

    def test_heatmap_invalid_input_rejected(self):
        cases = (pd.DataFrame(), pd.DataFrame([[1, 0]]), pd.DataFrame([[1, .1], [.2, 1]]),
                 pd.DataFrame([[1, np.inf], [np.inf, 1]]), pd.DataFrame([[1.01]]),
                 pd.DataFrame([[-.1]]), pd.DataFrame([[1]], index=['A'], columns=['B']))
        for frame in cases:
            with self.assertRaises(ValueError):
                ldheatmap.plot(frame)
        with self.assertRaises(ValueError):
            ldheatmap.plot(pd.DataFrame([[1]]), plot_diag=False)
        index = pd.MultiIndex.from_tuples([('chr1', 10), ('chr2', 20)], names=['chrom', 'pos'])
        with self.assertRaisesRegex(ValueError, 'one chromosome'):
            ldheatmap.plot(pd.DataFrame([[1, 0], [0, 1]], index=index, columns=index), plot_snp=True)


if __name__ == '__main__':
    unittest.main()
