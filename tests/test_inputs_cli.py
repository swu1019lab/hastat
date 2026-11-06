"""Input identifiers and real CLI integration; execute inside a PBS allocation."""
from functools import partial
from hastat.utils import population
import importlib.metadata
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd
import pysam

from hastat import __version__
from hastat.func import stat, view, plot
from hastat.utils.variant import VariantData
from hastat.dataio.table import read_assignments, read_sample_csv
from hastat.utils.association import HaplotypeAssociation, PAIR_COLUMNS, ASSOCIATION_VERSION
from hastat.viz.network import FancyNetwork
from test_population import write_vcf


class InputCLIValidation(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def write(self, name, text):
        path = self.root / name
        path.write_text(text, encoding='utf-8')
        return path

    def cli(self, *args):
        command = ([os.environ['HASTAT_CLI']] if os.environ.get('HASTAT_CLI')
                   else [sys.executable, '-m', 'hastat.hastat'])
        return subprocess.run(command + list(map(str, args)), capture_output=True, text=True)

    def test_missing_subcommand_exits_two_without_traceback(self):
        for args in ((), ('plot',)):
            with self.subTest(args=args):
                result = self.cli(*args)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertIn('usage:', result.stderr)
                self.assertNotIn('Traceback', result.stderr)

    def test_help_and_version(self):
        for args in (('--help',), ('plot', '--help'), ('view', '--help'),
                     ('stat', '--help'), ('network', '--help'), ('plot', 'network', '--help')):
            result = self.cli(*args)
            self.assertEqual(result.returncode, 0, result.stderr)
        help_text = self.cli('view', '--help').stdout
        self.assertIn('allow all heterozygosity rates', help_text)
        self.assertNotIn('only heterozygous sites', help_text)
        result = self.cli('--version')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip().split()[-1], __version__)
        self.assertEqual(ASSOCIATION_VERSION, __version__)
        if os.environ.get('HASTAT_CLI'):
            self.assertEqual(importlib.metadata.version('hastat'), __version__)

    def test_plot_direct_call_requires_type(self):
        with self.assertRaisesRegex(ValueError, 'No plot type'):
            plot.run(SimpleNamespace(plot_type=None))

    def test_two_column_assignments_keep_literal_ids_and_labels(self):
        path = self.write('groups.csv', '\ufeffsamples,groups\n001,01\n1,NA\n NA ,01\n')
        result = read_assignments(path, available_samples=['001', '1', 'NA'])
        self.assertEqual(result.samples.tolist(), ['001', '1', 'NA'])
        self.assertEqual(result.groups.tolist(), ['01', 'NA', '01'])

    def test_assignments_reject_malformed_or_duplicate_inputs(self):
        invalid = ('samples,groups\n', 'samples,groups\ns,A\n s ,B\n',
                   'samples,groups\n,A\n', 'samples,groups\ns,\n',
                   'samples,groups\nunknown,A\n', 'samples,groups,x\ns,A,1\n',
                   'samples,samples\ns,A\n')
        for text in invalid:
            with self.subTest(text=text), self.assertRaises(ValueError):
                read_assignments(self.write('bad.csv', text), available_samples=['s'])

    def stat_files(self):
        hap = self.write('hap.csv', 'samples,haplotypes\n001,A\n1,A\nNA,A\n002,B\n2,B\ns6,B\ns7,A\n')
        phe = self.write('phe.csv', 'samples,trait\n001,1\n1,2\nNA,4\n002,5\n2,7\ns6,8\ns7,NA\n')
        return hap, phe

    def test_stat_literal_na_id_and_missing_measurement(self):
        hap, phe = self.stat_files()
        output = self.root / 'result.csv'
        before = set(self.root.iterdir())
        result = self.cli('stat', '-g', hap, '-p', phe, '-s', 2, '-m', 'WelchBonferroni', '-o', output)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(set(self.root.iterdir()) - before, {output})
        frame = pd.read_csv(output)
        self.assertEqual(frame.columns.tolist(), PAIR_COLUMNS)
        self.assertEqual(frame[['count1', 'count2']].iloc[0].tolist(), [3, 3])
        reference = HaplotypeAssociation(pd.DataFrame({'samples': ['a', 'b', 'c', 'd', 'e', 'f'],
                                              'trait': [1, 2, 4, 5, 7, 8]}),
            pd.DataFrame({'samples': ['a', 'b', 'c', 'd', 'e', 'f'], 'haplotypes': ['A']*3 + ['B']*3}),
            min_hap_size=2, method='WelchBonferroni').comparisons
        np.testing.assert_allclose(frame[['mean1', 'mean2', 'p-adj']], reference[['mean1', 'mean2', 'p-adj']])

    def test_stat_duplicate_headers_and_samples_fail_before_output(self):
        hap, phe = self.stat_files()
        for text in ('samples,trait,trait\n001,1,2\n', 'samples,trait\n001,1\n 001 ,2\n'):
            phe.write_text(text)
            result = self.cli('stat', '-g', hap, '-p', phe, '-s', 2, '-o', self.root / 'invalid.csv')
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((self.root / 'invalid.csv').exists())

    def population_files(self):
        rows = [(1, 'A', 'G', ['0/0', '0/1', '1/1', '0/0']),
                (2, 'A', 'G', ['0/1', '0/0', '1/1', '0/1']),
                (3, 'A', 'G', ['1/1', '0/0', '0/1', '0/0'])]
        path = write_vcf(self.root / 'input.vcf', rows, samples=('001', '1', '002', '2'))
        groups = self.write('pop.csv', 'samples,groups\n001,NA\n1,NA\n002,01\n2,01\n')
        return path, groups

    def test_pi_and_d_numeric_ids_and_na_population_match_direct_api(self):
        path, groups = self.population_files()
        with pysam.VariantFile(path) as handle:
            gv = VariantData(handle)
            for kind, method, columns in (
                    ('pi', partial(population.calculate_pi, gv), ['chrom', 'start', 'end', 'n_bases', 'counts', 'pi', 'groups']),
                    ('tajima_d', partial(population.calculate_tajima_d, gv), ['chrom', 'start', 'end', 'counts', 'tajima_d', 'groups'])):
                result = self.cli('view', '-v', path, '-r', 'chr1:1-3', '-g', groups, '-t', kind,
                                  '--size', 3, '--step', 3, '-o', self.root / kind)
                self.assertEqual(result.returncode, 0, result.stderr)
                frame = pd.read_csv(self.root / f'{kind}.{kind}.csv', dtype={'groups': str}, keep_default_na=False)
                self.assertEqual(frame.columns.tolist(), columns)
                self.assertEqual(frame.groups.tolist(), ['NA', '01'])
                expected = [float(method('chr1', 1, 3, size=3, step=3, sample_list=ids)[0, -1])
                            for ids in (['001', '1'], ['002', '2'])]
                np.testing.assert_allclose(frame[kind].astype(float), expected)

    def test_population_invalid_assignments_fail_before_output(self):
        path, groups = self.population_files()
        for kind in ('pi', 'tajima_d', 'fst', 'freq'):
            for text in ('samples,groups\n001,A\n 001 ,B\n', 'samples,groups\nmissing,A\n1,B\n',
                         'samples,groups\n', 'samples,groups\n001,\n1,B\n'):
                groups.write_text(text)
                args = SimpleNamespace(vcf=[path], region='chr1:1-3', upstream=0, downstream=0,
                    type=kind, out=str(self.root / 'bad'), het=None, maf=0, max_missing=1,
                    group=str(groups), size=3, step=3)
                with self.subTest(kind=kind, text=text), self.assertRaises(ValueError):
                    view.run(args)
                self.assertFalse((self.root / f'bad.{kind}.csv').exists())

    def network_config(self):
        graph = self.write('net.csv', 'source,target,weight,type\nHap1,Hap2,1,edge\n')
        hap = self.write('nhap.csv', 'samples,haplotypes\n001,Hap1\n1,Hap2\nNA,Hap1\n')
        groups = self.write('ngroups.csv', 'samples,pop,region\n001,NA,N\n1,P,S\nNA,P,N\n')
        return {'data': {'file': str(graph), 'sample_hap': str(hap), 'sample_group': str(groups)}, 'plot': {}}

    def test_network_literal_ids_exact_composition_and_second_group(self):
        config = self.network_config()
        net = FancyNetwork(config)
        self.assertEqual(net.node_sizes, {'Hap1': 2, 'Hap2': 1})
        self.assertEqual(net.hap_composition, {'Hap1': {'NA': .5, 'P': .5}, 'Hap2': {'P': 1}})
        config['data']['group_index'] = 2
        net = FancyNetwork(config)
        self.assertEqual(net.hap_composition, {'Hap1': {'N': 1}, 'Hap2': {'S': 1}})

    def test_network_duplicate_ids_in_either_table_rejected(self):
        for table in ('sample_hap', 'sample_group'):
            config = self.network_config()
            path = Path(config['data'][table])
            path.write_text(path.read_text() + path.read_text().splitlines()[1] + '\n')
            with self.subTest(table=table), self.assertRaisesRegex(ValueError, 'duplicate'):
                FancyNetwork(config)

    def test_network_invalid_group_index_or_disjoint_samples_rejected(self):
        config = self.network_config()
        for index in (0, -1, 3, True, 1.5):
            config['data']['group_index'] = index
            with self.subTest(index=index), self.assertRaisesRegex(ValueError, 'group_index'):
                FancyNetwork(config)
        config['data']['group_index'] = 1
        Path(config['data']['sample_group']).write_text('samples,pop\nother,P\n')
        with self.assertRaisesRegex(ValueError, 'No overlapping'):
            FancyNetwork(config)

    def test_network_partial_join_warns_and_counts_matched_samples(self):
        config = self.network_config()
        Path(config['data']['sample_group']).write_text('samples,pop\n001,NA\n1,P\nother,P\n')
        with self.assertLogs('hastat', level='WARNING'):
            net = FancyNetwork(config)
        self.assertEqual(net.node_sizes, {'Hap1': 1, 'Hap2': 1})

    def test_network_without_optional_metadata(self):
        config = self.network_config()
        del config['data']['sample_group']
        self.assertEqual(FancyNetwork(config).node_sizes, {'Hap1': 2, 'Hap2': 1})
        del config['data']['sample_hap']
        self.assertEqual(FancyNetwork(config).node_sizes, {'Hap1': 1, 'Hap2': 1})

    def test_network_filters_only_missing_selected_metadata(self):
        config = self.network_config()
        Path(config['data']['sample_group']).write_text('samples,pop,other\n001,A,\n1,B,\nNA,,C\n')
        with self.assertLogs('hastat', level='WARNING'):
            net = FancyNetwork(config)
        self.assertEqual(net.node_sizes, {'Hap1':1,'Hap2':1})
        config['data']['group_index'] = 2
        with self.assertLogs('hastat', level='WARNING'):
            self.assertEqual(FancyNetwork(config).node_sizes, {'Hap1':1,'Hap2':0})

    def test_network_real_cli_dispatch_and_invalid_input_exit(self):
        config = self.network_config()['data']
        output = self.root / 'network.pdf'
        args = ['plot', 'network', '--file', config['file'], '--sample_hap', config['sample_hap'],
                '--sample_group', config['sample_group'], '-o', output]
        result = self.cli(*args)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertGreater(output.stat().st_size, 1000)
        output.unlink()
        path = Path(config['sample_hap'])
        path.write_text(path.read_text() + '001,Hap2\n')
        result = self.cli(*args)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
