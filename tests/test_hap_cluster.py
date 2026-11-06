"""Prefix discovery, category validation and independent clustering fixtures."""
import itertools
from pathlib import Path
import tempfile
import unittest

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from hastat.utils.haplotype import cluster_haplotypes
from hastat.viz.gene import FancyGene
from test_hap_counts import example_inputs, example_config, labels_for


class ClusterMathValidation(unittest.TestCase):
    def fixture(self):
        # Hand-calculated Hamming: AB=1/4, AC=1, BC=3/4.
        return pd.DataFrame({'A':[0,0,0,0], 'B':[0,0,0,2], 'C':[2,2,2,2]})

    def test_average_matches_hand_calculated_merges(self):
        names, z = cluster_haplotypes(self.fixture(), ['A','B','C'])
        np.testing.assert_allclose(z, [[0,1,.25,2],[2,3,.875,3]])
        self.assertEqual(names, ['C','A','B'])

    def test_other_supported_linkage_heights(self):
        for method, height in [('single',.75),('complete',1.),('weighted',.875)]:
            with self.subTest(method=method):
                _, z = cluster_haplotypes(self.fixture(), ['A','B','C'],method)
                np.testing.assert_allclose(z[:,2], [.25,height])

    def test_input_order_and_tied_distances_are_stable(self):
        table = pd.DataFrame({'A':[0,0,2], 'B':[0,2,0], 'C':[2,0,0], 'D':[2,2,2]})
        expected_names, expected_z = cluster_haplotypes(table,list(table))
        for names in itertools.permutations(table.columns):
            actual_names, actual_z = cluster_haplotypes(table[list(names)].iloc[::-1],list(names))
            self.assertEqual(actual_names, expected_names)
            np.testing.assert_array_equal(actual_z, expected_z)

    def test_missing_loci_excluded_for_all_pairs(self):
        table = self.fixture()
        table.loc[4] = [0,2,-1]
        with self.assertLogs('hastat',level='WARNING'):
            _, z = cluster_haplotypes(table,list(table))
        np.testing.assert_allclose(z[:,2],[.25,.875])

    def test_no_common_loci_and_invalid_genotypes_fail(self):
        with self.assertRaisesRegex(ValueError,'No complete loci'):
            cluster_haplotypes(pd.DataFrame({'A':[0,-1],'B':[-1,2]}),['A','B'])
        for value in [3,np.nan,np.inf]:
            with self.subTest(value=value), self.assertRaisesRegex(ValueError,'genotypes'):
                cluster_haplotypes(pd.DataFrame({'A':[0],'B':[value]}),['A','B'])

    def test_unsupported_methods_fail_instead_of_omitting_tree(self):
        for method in ['ward','centroid','median','typo']:
            with self.subTest(method=method), self.assertRaisesRegex(ValueError,'Euclidean'):
                cluster_haplotypes(self.fixture(),['A','B','C'],method)

    def test_two_one_and_identical_profiles(self):
        _, z = cluster_haplotypes(self.fixture(),['B','A'])
        np.testing.assert_allclose(z,[[0,1,.25,2]])
        self.assertEqual(cluster_haplotypes(self.fixture(),['A']),(['A'],None))
        names, z = cluster_haplotypes(pd.DataFrame({'A':[0,0],'B':[0,0],'C':[0,0]}),['C','B','A'])
        self.assertEqual(set(names),{'A','B','C'})
        np.testing.assert_array_equal(z[:,2],[0,0])


class PrefixPlotValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.addCleanup(plt.close,'all')
        self.root = Path(self.tmp.name)
        example_inputs(self.root)

    def draw(self, extra=''):
        obj = FancyGene(example_config(self.root,extra))
        obj.plot()
        return obj

    def test_absolute_prefix_and_explicit_table_path(self):
        for prefix in [str(self.root/'gene'),str(self.root/'gene.table.csv')]:
            config = example_config(self.root)
            path = self.root/'gene.toml'
            path.write_text(path.read_text().replace('hap_file="gene"',f'hap_file="{prefix}"'))
            obj = FancyGene(config)
            obj.plot()
            self.assertEqual(set(labels_for(obj)),{'Hap1 (n=25)','Hap2 (n=12)','Hap3 (n=3)'})

    def test_absent_table_error_and_exact_prefix_pairing(self):
        (self.root/'gene.table.csv').rename(self.root/'other.table.csv')
        with self.assertRaises(FileNotFoundError):
            self.draw()
        (self.root/'other.table.csv').rename(self.root/'gene.table.csv')
        (self.root/'gene.group.csv').rename(self.root/'other.group.csv')
        self.assertEqual(labels_for(self.draw()),['Hap1','Hap2','Hap3'])

    def test_present_malformed_group_rejected_even_if_counts_hidden(self):
        (self.root/'gene.group.csv').write_text('samples\ns1\n')
        for extra in ['', 'hap_count_show=false\n']:
            with self.subTest(extra=extra), self.assertRaises(ValueError):
                self.draw(extra)

    def test_unknown_assignments_rejected(self):
        path = self.root/'gene.group.csv'
        path.write_text(path.read_text()+'other,Hap99\n')
        with self.assertRaisesRegex(ValueError,'Hap99'):
            self.draw()

    def test_unknown_keys_warn_and_plot_without_interpreting_their_values(self):
        keys = ['xy', 'width', 'hap_group_file', 'misspelled_parameter']
        with self.assertLogs('hastat', level='WARNING') as captured:
            obj = self.draw(''.join(f'{key}="ignored"\n' for key in keys))
        messages = [record.getMessage() for record in captured.records]
        for key in keys:
            self.assertIn(f"Unknown gene TOML parameter '{key}' in [default]", messages)
        self.assertEqual(set(labels_for(obj)), {'Hap1 (n=25)', 'Hap2 (n=12)', 'Hap3 (n=3)'})

    def test_invalid_color_assignments_fail(self):
        for text in ['Hap1,A\nHap1,B\n','Hap1,A\nHap1,A\n','Hap1,\n','Hap1,A,extra\n']:
            (self.root/'colors.csv').write_text(text)
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.draw('hap_color_file="colors.csv"\n')

    def test_table_identifiers_values_and_duplicate_headers(self):
        obj = FancyGene({'plot':{}})
        path = self.root/'bad.csv'
        for text in ['chrom,pos,ref,alt,Hap1,Hap1\n1,5,A,G,0,2\n',
                     'chrom,pos,ref,alt,Hap1\n1,0,A,G,0\n',
                     'chrom,pos,ref,alt,Hap1\n1,5,A,G,3\n',
                     'chrom,pos,ref,alt,Hap1\n1,5,A,G,abc\n',
                     'chrom,pos,ref,alt,Hap1\n1,5,A,G,0\n1,5,A,G,2\n']:
            path.write_text(text)
            with self.subTest(text=text), self.assertRaises(ValueError):
                obj._load_haplotype_data(path)
        path.write_text('chrom,pos,ref,alt,Hap1,Hap2\nNA,5,A,G,NA,\n')
        frame = obj._load_haplotype_data(path)
        self.assertEqual(frame.chrom.iloc[0],'NA')
        self.assertEqual(frame[['Hap1','Hap2']].iloc[0].tolist(),[-1,-1])

    def test_two_haplotypes_render_tree_with_or_without_counts(self):
        for counts in [True,False]:
            if not counts:
                (self.root/'gene.group.csv').unlink()
            obj = self.draw('hap_cols=["Hap1","Hap2"]\nhap_cluster_show=true\n')
            labels,bars,tree = obj._count_label_layout[0]
            self.assertEqual(len(labels),2)
            self.assertEqual(len(tree),1)
            obj.fig.canvas.draw()
            renderer = obj.fig.canvas.get_renderer()
            self.assertLess(max(x.get_window_extent(renderer).x1 for x in tree),
                            min(x.get_window_extent(renderer).x0 for x in labels))

    def test_bad_cluster_method_and_no_common_loci_abort_plot(self):
        with self.assertRaisesRegex(ValueError,'Euclidean'):
            self.draw('hap_cluster_show=true\nhap_cluster_linkage="ward"\n')
        path = self.root/'gene.table.csv'
        frame = pd.read_csv(path)
        frame['Hap1'] = -1
        frame.to_csv(path,index=False)
        with self.assertRaisesRegex(ValueError,'No complete loci'):
            self.draw('hap_cluster_show=true\n')


if __name__ == '__main__':
    unittest.main(verbosity=2)
