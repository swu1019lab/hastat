"""Network algorithm and export integration checks; run in a PBS allocation."""
from hastat.dataio.network import save_network
import argparse
import ast
from importlib.util import find_spec
import contextlib
import io
import itertools
import logging
from pathlib import Path
import random
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import networkx as nx
import numpy as np
import pandas as pd

from hastat.func import network as network_command
from hastat.log.logger import logger
from hastat.utils.network import HaplotypeNetwork
from hastat.dataio.network import read_network
from hastat.viz.network import FancyNetwork

ROOT = Path(__file__).resolve().parents[1]
SQUARE = pd.DataFrame({'A': [0, 0], 'B': [0, 1], 'C': [1, 0], 'D': [1, 1]})
CHAIN = pd.DataFrame({'A': [0, 0, 0], 'B': [0, 0, 1], 'C': [0, 1, 1], 'D': [1, 1, 1]})


def edge_set(graph):
    return {frozenset((u, v)) for u, v in graph.edges()}


def weighted_edges(graph):
    return {frozenset((u, v)): d['weight'] for u, v, d in graph.edges(data=True)}


def make_network(frame, method='MSN'):
    net = HaplotypeNetwork()
    net.set_haplotypes(frame)
    net.create_network(method)
    return net


def brute_force_tree_union(nodes, edges):
    """Enumerate edge subsets and check connectivity by traversal, without MST code."""
    best, union, optimum_count = float('inf'), set(), 0
    for chosen in itertools.combinations(edges, len(nodes) - 1):
        adjacent = {node: set() for node in nodes}
        for u, v, _ in chosen:
            adjacent[u].add(v)
            adjacent[v].add(u)
        reached, pending = set(), [nodes[0]]
        while pending:
            node = pending.pop()
            if node not in reached:
                reached.add(node)
                pending.extend(adjacent[node] - reached)
        if len(reached) != len(nodes):
            continue
        total = sum(w for _, _, w in chosen)
        selected = {frozenset((u, v)) for u, v, _ in chosen}
        if total < best:
            best, union, optimum_count = total, selected, 1
        elif total == best:
            union.update(selected)
            optimum_count += 1
    return best, union, optimum_count


def write_table(path, frame):
    metadata = pd.DataFrame({'chrom': ['chr1'] * len(frame), 'pos': range(1, len(frame) + 1),
                             'ref': ['A'] * len(frame), 'alt': ['G'] * len(frame)})
    pd.concat([metadata, frame.reset_index(drop=True)], axis=1).to_csv(path, index=False)


class NetworkValidation(unittest.TestCase):
    def assert_oracle(self, nodes, edges):
        minimum, union, _ = brute_force_tree_union(nodes, edges)
        net = HaplotypeNetwork()
        net.nodes, net.edge_data = nodes, edges
        net.create_network('MSN')
        self.assertEqual(edge_set(net.network), union)
        net.create_network('MST')
        self.assertTrue(nx.is_tree(net.network))
        self.assertEqual(sum(w for _, _, w in net.network.edges(data='weight')), minimum)

    def test_exhaustive_k4_weights(self):
        nodes = list('ABCD')
        pairs = list(itertools.combinations(nodes, 2))
        with patch.object(logger, 'info'):
            for weights in itertools.product(range(3), repeat=len(pairs)):
                self.assert_oracle(nodes, [(u, v, w) for (u, v), w in zip(pairs, weights)])

    def test_seeded_k5_weights(self):
        rng, nodes = random.Random(20261031), list('ABCDE')
        with patch.object(logger, 'info'):
            for _ in range(30):
                self.assert_oracle(nodes, [(u, v, rng.randrange(5)) for u, v in itertools.combinations(nodes, 2)])

    def test_square_retains_only_tied_sides(self):
        net = make_network(SQUARE)
        self.assertEqual(edge_set(net.network), {frozenset(e) for e in ('AB', 'AC', 'BD', 'CD')})
        self.assertTrue(all(w == 1 for _, _, w in net.network.edges(data='weight')))
        self.assertEqual(make_network(SQUARE, 'MST').network.number_of_edges(), 3)

    def test_chain_has_no_extra_edges(self):
        self.assertEqual(edge_set(make_network(CHAIN).network), {frozenset(e) for e in ('AB', 'BC', 'CD')})

    def test_equal_distances_can_legitimately_be_complete(self):
        triangle = pd.DataFrame({'A': [0, 0, 1], 'B': [0, 1, 0], 'C': [1, 0, 0]})
        graph = make_network(triangle).network
        self.assertEqual(graph.number_of_edges(), 3)
        self.assertTrue(all(w == 2 for _, _, w in graph.edges(data='weight')))

    def test_zero_distances_keep_distinct_nodes(self):
        frame = pd.DataFrame({'A': [0, 1], 'B': [0, 1], 'C': [0, 1]})
        net = make_network(frame)
        self.assertEqual(set(net.network), set('ABC'))
        self.assertEqual(net.network.number_of_edges(), 3)
        self.assertEqual(sum(w for _, _, w in net.network.edges(data='weight')), 0)
        self.assertEqual(make_network(frame, 'MST').network.number_of_edges(), 2)

    def test_missing_loci_are_excluded_for_every_pair(self):
        frame = pd.DataFrame({'A': [0, -1, 0, 2], 'B': [1, 1, np.nan, 2], 'C': [1, 0, 1, 0]})
        self.assertEqual(HaplotypeNetwork().calculate_haplotype_distances(frame), [('A', 'B', 1), ('A', 'C', 2), ('B', 'C', 1)])
        # Even the two observed entries at an excluded locus must not affect distance.
        frame.loc[1, ['B', 'C']] = [2, 2]
        frame.loc[2, ['A', 'C']] = [2, 2]
        self.assertEqual(HaplotypeNetwork().calculate_haplotype_distances(frame), [('A', 'B', 1), ('A', 'C', 2), ('B', 'C', 1)])

    def test_state_distance_is_not_dosage_distance(self):
        frame = pd.DataFrame({'A': [0], 'B': [2]})
        self.assertEqual(HaplotypeNetwork().calculate_haplotype_distances(frame), [('A', 'B', 1)])

    def test_no_complete_loci_rejected(self):
        for frame in (pd.DataFrame({'A': []}), pd.DataFrame({'A': [-1]}),
                      pd.DataFrame({'A': [0, np.nan], 'B': [-1, 0]})):
            with self.subTest(frame=frame), self.assertRaisesRegex(ValueError, 'No loci'):
                make_network(frame)

    def test_invalid_states_and_names_rejected(self):
        for value in (3, -2, 0.5, float('inf'), 'wrong'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                make_network(pd.DataFrame({'A': [0], 'B': [value]}))
        for frame in (pd.DataFrame([[0, 1]], columns=['A', 'A']), pd.DataFrame({'': [0]}),
                      pd.DataFrame({'A\tB': [0]}), pd.DataFrame({None: [0]}), pd.DataFrame()):
            with self.subTest(columns=frame.columns), self.assertRaises(ValueError):
                make_network(frame)

    def test_column_and_locus_order_invariant(self):
        for method in ('MST', 'MSN'):
            reference = weighted_edges(make_network(SQUARE, method).network)
            for columns in itertools.permutations(SQUARE.columns):
                graph = make_network(SQUARE.loc[::-1, list(columns)], method).network
                self.assertEqual(weighted_edges(graph), reference)

    def test_disconnected_forest_and_isolate(self):
        net = HaplotypeNetwork()
        net.nodes = list('ABCDEF')
        net.edge_data = [('A', 'B', 1), ('B', 'C', 1), ('A', 'C', 1), ('D', 'E', 2)]
        net.create_network('MSN')
        self.assertEqual(net.network.number_of_edges(), 4)
        self.assertEqual(set(net.network), set('ABCDEF'))
        net.create_network('MST')
        self.assertEqual(net.network.number_of_edges(), 3)
        self.assertTrue(nx.is_forest(net.network))

    def test_invalid_direct_edges_rejected(self):
        for edges in ([('A', 'A', 0)], [('A', 'B', -1)], [('A', 'B', np.nan)],
                      [('A', 'B', 1), ('B', 'A', 2)]):
            net = HaplotypeNetwork()
            net.edge_data = edges
            with self.assertRaises(ValueError):
                net.create_network('MSN')

    def test_singleton_all_formats_and_plot_loader(self):
        with tempfile.TemporaryDirectory() as tmp:
            for method in ('MST', 'MSN'):
                net = make_network(pd.DataFrame({'Hap1': [0, 1]}), method)
                self.assertEqual(list(net.network), ['Hap1'])
                self.assertEqual(net.network.number_of_edges(), 0)
                for suffix in ('csv', 'txt', 'gexf', 'gml'):
                    with self.subTest(method=method, suffix=suffix):
                        path = Path(tmp) / f'{method}.{suffix}'
                        save_network(net.network, str(path))
                        graph = read_network(path)
                        self.assertEqual(list(graph), ['Hap1'])
                        self.assertEqual(graph.number_of_edges(), 0)
                        plot = FancyNetwork({'data': {'file': str(path)}, 'plot': {}})
                        self.assertEqual(plot.node_sizes, {'Hap1': 1})
                frame = pd.read_csv(Path(tmp) / f'{method}.csv', keep_default_na=False)
                self.assertEqual(list(frame.columns), ['source', 'target', 'weight', 'type'])
                self.assertEqual(frame.iloc[0].tolist(), ['Hap1', '', '', 'node'])

    def test_regular_round_trip_and_csv_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            net = make_network(SQUARE)
            for suffix in ('csv', 'txt', 'gexf', 'gml'):
                path = Path(tmp) / f'network.{suffix}'
                save_network(net.network, str(path))
                self.assertEqual(weighted_edges(read_network(path)), weighted_edges(net.network))
                self.assertEqual(set(read_network(path)), set(net.network))
            frame = pd.read_csv(Path(tmp) / 'network.csv')
            self.assertEqual(list(frame.columns), ['source', 'target', 'weight', 'type'])
            self.assertEqual(set(frame.type), {'edge'})
            self.assertFalse(list(Path(tmp).glob('*.summary.json')))

    def test_zero_weight_edges_and_isolate_round_trip(self):
        net = HaplotypeNetwork()
        net.nodes, net.edge_data = ['A', 'B', 'C'], [('A', 'B', 0)]
        net.create_network('MSN')
        with tempfile.TemporaryDirectory() as tmp:
            for suffix in ('csv', 'txt', 'gexf', 'gml'):
                path = Path(tmp) / f'network.{suffix}'
                save_network(net.network, str(path))
                self.assertEqual(set(read_network(path)), set('ABC'))
                self.assertEqual(weighted_edges(read_network(path)), {frozenset('AB'): 0})

    def test_literal_node_names_survive(self):
        with tempfile.TemporaryDirectory() as tmp:
            for label in ('NA', '001', '#node', 'A#B', '甲'):
                for frame in (pd.DataFrame({label: [0]}), pd.DataFrame({label: [0], 'other': [1]})):
                    net = make_network(frame)
                    for suffix in ('csv', 'txt', 'gexf', 'gml'):
                        path = Path(tmp) / f'network.{suffix}'
                        save_network(net.network, str(path))
                        self.assertEqual(set(read_network(path)), set(frame.columns))

    def test_command_callback_and_invalid_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp) / 'input.csv', Path(tmp) / 'result.csv'
            write_table(source, SQUARE)
            args = SimpleNamespace(hap_table=str(source), method='MSN', out=str(output))
            network_command.run(args)
            self.assertEqual(weighted_edges(read_network(output)), weighted_edges(make_network(SQUARE).network))
            self.assertEqual(set(p.name for p in Path(tmp).iterdir()), {'input.csv', 'result.csv'})
            for text in ('wrong,pos,ref,alt,A\nchr1,1,A,G,0\n',
                         'chrom,pos,ref,alt,A,A\nchr1,1,A,G,0,1\n',
                         'chrom,pos,ref,alt\nchr1,1,A,G\n'):
                source.write_text(text)
                with self.assertRaises(ValueError):
                    network_command.run(args)

    def test_cli_contract_and_actual_network_dispatch(self):
        tree = ast.parse(Path(find_spec('hastat.hastat').origin).read_text(encoding='utf-8'))
        main = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
        calls = []
        fake = SimpleNamespace(run=lambda args: calls.append(args.method))
        namespace = dict(argparse=argparse, logging=logging, logger=logging.getLogger('cli_contract'),
                         view=fake, stat=fake, plot=fake, gwas=fake, network=network_command)
        exec(compile(ast.Module(body=[main], type_ignores=[]), 'hastat/hastat.py', 'exec'), namespace)
        with tempfile.TemporaryDirectory() as tmp:
            source, output = Path(tmp) / 'input.csv', Path(tmp) / 'network.csv'
            write_table(source, SQUARE)
            for method, count in (('MST', 3), ('MSN', 4)):
                args = ['hastat', 'network', '-t', str(source), '-m', method, '-o', str(output)]
                with patch.object(sys, 'argv', args):
                    namespace['main']()
                self.assertEqual(read_network(output).number_of_edges(), count)
            for removed in ('-f', '--threshold_factor'):
                with patch.object(sys, 'argv', args + [removed, '1.2']), contextlib.redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit) as exc:
                        namespace['main']()
                    self.assertEqual(exc.exception.code, 2)
            for method in ('TukeyHSD', 'WelchBonferroni'):
                with patch.object(sys, 'argv', ['hastat', 'stat', '-g', 'groups', '-p', 'pheno', '-m', method]):
                    namespace['main']()
            self.assertEqual(calls, ['TukeyHSD', 'WelchBonferroni'])

    def test_invalid_api_use_rejected(self):
        with self.assertRaises(ValueError):
            save_network(HaplotypeNetwork().network, 'unused.csv')
        with self.assertRaises(ValueError):
            HaplotypeNetwork().create_network()
        net = make_network(SQUARE)
        with self.assertRaises(ValueError):
            net.create_network('unknown')
        with self.assertRaises(TypeError):
            net.create_network('MSN', threshold_factor=1.2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
