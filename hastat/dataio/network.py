"""Read hastat network formats, including explicitly stored isolated nodes."""
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from hastat.log.logger import logger


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


def save_network(graph, file_path: str, format: str = 'auto'):
    """
    Save the network to a file in various formats

    :param file_path: The path to save the network
    :param format: File format ('gexf', 'gml', 'edgelist', 'csv', 'auto').
                  If 'auto', format is determined by file extension.
    """
    from pathlib import Path

    if graph is None:
        raise ValueError("No network created yet. Please call create_network() first.")

    # Determine format from file extension if auto
    if format == 'auto':
        ext = Path(file_path).suffix.lower()
        format_map = {
            '.gexf': 'gexf',
            '.gml': 'gml',
            '.txt': 'edgelist',
            '.edgelist': 'edgelist',
            '.csv': 'csv'
        }
        format = format_map.get(ext, 'csv')

    logger.info(f"Saving network in {format.upper()} format to: {file_path}")

    # Save in requested format
    if format == 'gexf':
        # Add node attributes for better visualization
        for node in graph.nodes():
            graph.nodes[node]['label'] = node
            graph.nodes[node]['id'] = node

        # Add edge attributes
        for u, v, data in graph.edges(data=True):
            data['label'] = f"weight: {data['weight']}"
            data['type'] = 'undirected'

        nx.write_gexf(graph, file_path)

    elif format == 'gml':
        # GML format - ensure all attributes are valid
        network_copy = graph.copy()

        # Add node attributes
        for node in network_copy.nodes():
            network_copy.nodes[node]['label'] = str(node)
            network_copy.nodes[node]['id'] = str(node)

        # Ensure edge weights are numeric
        for u, v, data in network_copy.edges(data=True):
            if 'weight' in data:
                data['weight'] = float(data['weight'])

        nx.write_gml(network_copy, file_path)

    elif format == 'edgelist':
        # Simple edge list format
        import json
        with open(file_path, 'w', encoding='utf-8') as f:
            for node in nx.isolates(graph):
                # A comment keeps normal edge-list readers compatible;
                # hastat's reader restores the otherwise unrepresentable node.
                f.write('#node\t' + json.dumps(str(node), ensure_ascii=False) + '\n')
            for u, v, data in graph.edges(data=True):
                weight = data.get('weight', 1)
                f.write(f"{u}\t{v}\t{weight}\n")

    elif format == 'csv':
        # Original CSV format with additional information
        results = []

        # Save edge information from the constructed network
        for u, v, data in graph.edges(data=True):
            results.append({
                'source': u,
                'target': v,
                'weight': data.get('weight', 0),
                'type': 'edge'
            })

        for node in nx.isolates(graph):
            results.append({'source': node, 'target': '', 'weight': '', 'type': 'node'})
        # Existing columns also represent isolated nodes, without self-loops.
        results_df = pd.DataFrame(results, columns=['source', 'target', 'weight', 'type'])
        results_df.to_csv(file_path, index=False)

    else:
        raise ValueError(f"Unsupported format: {format}. Supported formats: gexf, gml, edgelist, csv")

    logger.info(f"Network saved successfully with {graph.number_of_nodes()} nodes and {graph.number_of_edges()} edges")
