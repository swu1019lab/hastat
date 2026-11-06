"""Haplotype distances and minimum-spanning networks."""
import numpy as np
import pandas as pd
import networkx as nx
from hastat.log.logger import logger


class HaplotypeNetwork:
    """
    A class for gene haplotype network analysis
    """
    def __init__(self):
        self.edge_data = None
        self.network = None
        self.nodes = []

    @staticmethod
    def _node_key(node):
        return type(node).__name__, str(node)

    def calculate_haplotype_distances(self, df: pd.DataFrame):
        """
        Calculate the Hamming distance for all pairs of haplotypes

        :param df: DataFrame with haplotypes. each row is a locus. each column is a haplotype. values are 0, 1 or 2.
        :return: a list of tuples, each tuple contains the source, target and weight of an edge
        """
        if not isinstance(df, pd.DataFrame) or df.shape[1] == 0:
            raise ValueError("Haplotype table must contain at least one haplotype")
        if not df.columns.is_unique or df.columns.hasnans:
            raise ValueError("Haplotype names must be unique and non-missing")
        if any(not str(name).strip() or any(c in str(name) for c in '\t\r\n') for name in df.columns):
            raise ValueError("Haplotype names must be non-empty and contain no tab/newline")
        haplotypes = sorted(df.columns, key=self._node_key)
        try:
            values = df.loc[:, haplotypes].apply(pd.to_numeric, errors='raise').to_numpy(dtype=float, na_value=np.nan)
        except (TypeError, ValueError) as exc:
            raise ValueError("Haplotype states must be 0, 1, 2, -1 or missing") from exc
        if not np.all(np.isnan(values) | np.isin(values, [-1, 0, 1, 2])):
            raise ValueError("Haplotype states must be 0, 1, 2, -1 or missing")
        # All pairs must use the same loci, preserving an ordinary Hamming
        # metric. Neither -1 nor NaN is a biological genotype state.
        complete = np.all(np.isfinite(values) & (values != -1), axis=1)
        if not np.any(complete):
            raise ValueError("No loci with complete calls across all haplotypes")
        excluded = int(np.count_nonzero(~complete))
        if excluded:
            logger.warning(f"Network distances exclude {excluded} of {len(values)} loci "
                           "with missing calls across haplotypes")
        values = values[complete].T
        return [(haplotypes[i], haplotypes[j], int(np.count_nonzero(values[i] != values[j])))
                for i in range(len(haplotypes)) for j in range(i + 1, len(haplotypes))]

    def set_haplotypes(self, df: pd.DataFrame):
        """
        Set the haplotype states and calculate pairwise distances.

        :param df: DataFrame with haplotypes. each row is a locus. each column is a haplotype. values are 0, 1 or 2.
        """
        edges = self.calculate_haplotype_distances(df)
        self.nodes = sorted(df.columns, key=self._node_key)
        self.edge_data = edges
        self.network = None

    def _distance_graph(self):
        """Build a graph in canonical order, including nodes with no edges."""
        from numbers import Real
        if self.edge_data is None:
            raise ValueError("No edge data available. Call set_haplotypes() first.")
        edges, nodes = {}, set(self.nodes)
        for u, v, weight in self.edge_data:
            if u == v:
                raise ValueError("Distance edges must connect distinct nodes")
            if not isinstance(weight, Real) or not np.isfinite(weight) or weight < 0:
                raise ValueError("Distance weights must be finite non-negative numbers")
            u, v = sorted((u, v), key=self._node_key)
            if (u, v) in edges and edges[u, v] != weight:
                raise ValueError("Conflicting duplicate distance edges")
            edges[u, v] = weight
            nodes.update((u, v))
        graph = nx.Graph()
        graph.add_nodes_from(sorted(nodes, key=self._node_key))
        ordered = sorted(edges.items(), key=lambda item: (item[1], self._node_key(item[0][0]), self._node_key(item[0][1])))
        graph.add_weighted_edges_from((u, v, weight) for (u, v), weight in ordered)
        return graph

    def minimum_spanning_tree(self):
        """
        Construct the minimum spanning tree from distance matrix
        """
        return nx.minimum_spanning_tree(self._distance_graph(), algorithm='kruskal', weight='weight')

    def minimum_spanning_network(self):
        """
        Return the union of all minimum spanning forests, retaining tied edges.

        Process equal-weight edges as a batch. An edge is eligible exactly
        when its endpoints are disconnected using strictly lighter edges.
        Union components only after every edge in the batch has been tested.
        """
        from itertools import groupby
        graph = self._distance_graph()
        network = nx.Graph()
        network.add_nodes_from(graph.nodes())
        components = nx.utils.UnionFind(graph.nodes())
        edges = sorted(graph.edges(data=True), key=lambda edge: edge[2]['weight'])
        for _, group in groupby(edges, key=lambda edge: edge[2]['weight']):
            batch = list(group)
            for u, v, data in batch:
                if components[u] != components[v]:
                    network.add_edge(u, v, weight=data['weight'])
            for u, v, _ in batch:
                components.union(u, v)
        logger.info(f"MSN contains {network.number_of_nodes()} nodes and "
                    f"{network.number_of_edges()} edges, including equal-weight alternatives")
        return network

    def create_network(self, method: str = "MST"):
        """
        Create the network from distance matrix using different methods

        :param method: Method to construct the network. Available methods:
                      - "MST": Minimum Spanning Tree
                      - "MSN": Minimum Spanning Network
        """
        if self.edge_data is None:
            raise ValueError("No edge data available. Call set_haplotypes() first.")

        logger.info(f"Creating network using method: {method}")

        if method.upper() == "MST":
            self.network = self.minimum_spanning_tree()
        elif method.upper() == "MSN":
            self.network = self.minimum_spanning_network()
        else:
            available_methods = ["MST", "MSN"]
            raise ValueError(f"Invalid method: {method}. Available methods: {available_methods}")
