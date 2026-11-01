"""Read hastat network formats, including explicitly stored isolated nodes."""
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


def read_network(file_path):
    suffix = Path(file_path).suffix.lower()
    if suffix == '.gexf':
        return nx.read_gexf(file_path)
    if suffix == '.gml':
        return nx.read_gml(file_path)
    if suffix in ('.txt', '.edgelist'):
        edge_lines, isolated_nodes = [], []
        with open(file_path, encoding='utf-8') as handle:
            for line in handle:
                if line.startswith('#node\t') and len(line.rstrip('\r\n').split('\t')) == 2:
                    node = json.loads(line.split('\t', 1)[1])
                    if not isinstance(node, str) or not node:
                        raise ValueError('Invalid isolated node in network edge list')
                    isolated_nodes.append(node)
                else:
                    edge_lines.append(line)
        graph = nx.parse_edgelist(edge_lines, data=(('weight', float),), delimiter='\t', comments=None)
        graph.add_nodes_from(isolated_nodes)
        return graph
    data = pd.read_csv(file_path, dtype={'source': str, 'target': str, 'type': str}, keep_default_na=False)
    if list(data.columns) != ['source', 'target', 'weight', 'type']:
        raise ValueError('Network CSV must contain source,target,weight,type')
    graph = nx.Graph()
    for row in data.itertuples(index=False):
        if not row.source:
            raise ValueError('Network CSV contains an empty source node')
        if row.type == 'node':
            if row.target != '' or row.weight != '':
                raise ValueError('Isolated-node CSV rows must have empty target and weight')
            graph.add_node(row.source)
        elif row.type == 'edge':
            if not row.target:
                raise ValueError('Network edge has an empty target node')
            weight = float(row.weight)
            if not np.isfinite(weight) or weight < 0:
                raise ValueError('Network edge weight must be finite and non-negative')
            graph.add_edge(row.source, row.target, weight=weight)
        else:
            raise ValueError(f'Unknown network CSV row type: {row.type}')
    return graph
