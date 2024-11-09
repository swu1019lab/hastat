# -*- coding: utf-8 -*-
# @Time    : 2024/10/8 16:01
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : network.py

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib as mpl
from matplotlib import axes, patches, path


class HapNetwork(object):
    def __init__(self):
        self.edge_data = None
        self.node_data = None

    def add_edge_from_hap(self, df: pd.DataFrame, subset: list = None):
        # a dataframe contained haplotype group and genotype data
        if not isinstance(df, pd.DataFrame):
            raise ValueError("df is not a DataFrame")
        if subset is not None:
            df = df[df.iloc[:, 0].isin(subset)]
        name = df.iloc[:, 0].values
        geno = df.iloc[:, 1:].values
        # check the name whether is unique
        if len(name) != len(set(name)):
            raise ValueError("haplotype name is not unique")
        # calculate the distance between each pair of haplotypes
        dist_edges = []
        for i in range(len(name)):
            for j in range(i + 1, len(name)):
                dist = np.sum(geno[i] != geno[j])
                dist_edges.append((name[i], name[j], dist))
        self.edge_data = dist_edges

    def add_edge(self, df: pd.DataFrame, source: str = 'source', target: str = 'target', weight: str = 'weight'):
        """
        Add a dataframe contained edge data

        :param df: a dataframe contained edge data, the format should be:
            1. two columns
            2. 1st column is source nodes
            3. 2nd column is target nodes
            4. 3rd column is weight of edge
            5. other columns are edge attributes
        :param source: the name of source column
        :param target: the name of target column
        :param weight: the name of weight column
        :return:
        """
        if not isinstance(df, pd.DataFrame):
            raise ValueError("df is not a DataFrame")
        if np.intersect1d(df.columns.tolist(), [source, target, weight]).size < 3:
            raise ValueError("df is not a correct DataFrame")

        self.edge_data = df

    def add_node(self, df: pd.DataFrame, node: str = 'node', node_pie_data: list = None):
        """
        Add a dataframe contained node data

        :param df: a dataframe contained node data, the format should be:
            1. two columns at least
            2. 1st column is node names
            3. other columns are node attributes
        :param node: the name of node column
        :param node_pie_data: the name of node pie data columns
        :return:
        """
        if not isinstance(df, pd.DataFrame):
            raise ValueError("df is not a DataFrame")
        if node not in df.columns:
            raise ValueError("node is not in the DataFrame")
        if node_pie_data is None:
            raise ValueError("node pie data is None")
        self.node_data = df.loc[:, [node] + node_pie_data]
        self.node_data.columns = ['node'] + node_pie_data

    def create_graph(self):
        if isinstance(self.edge_data, pd.DataFrame):
            G = nx.from_pandas_edgelist(self.edge_data, edge_attr=True, create_using=nx.Graph())
        elif isinstance(self.edge_data, list):
            G = nx.Graph()
            G.add_weighted_edges_from(self.edge_data)
        else:
            raise ValueError("edge_data is not a DataFrame or list")
        return G

    def plot(self, ax: axes.Axes,
             scale_node_size: float = 1,
             node_font_size: float = 12, node_pie_colors: list = None,
             edge_color: str = 'k', k: float = None, scale: float = 1,
             seed: int = None) -> None:
        """
        Plot network for haplotype data

        :param ax: axes object
        :param scale_node_size: scale factor for node size
        :param node_font_size: font size for node
        :param node_pie_colors: colors for each pie chart of node
        :param edge_color: color for edge
        :param k: Optimal distance between nodes, default is 1/sqrt(n) where n is the number of nodes
        :param scale: Scale factor for positions
        :param seed: Seed for random layout
        :return:
        """
        G = self.create_graph()
        # Find the minimum spanning tree
        T = nx.minimum_spanning_tree(G)
        # plot the node labels
        pos = nx.spring_layout(T, k=k, scale=scale, seed=seed)
        nx.draw_networkx_labels(T, pos, font_size=node_font_size, font_family='arial')
        # plot the nodes with scatter-pie symbol, node attribute will be shown in pie chart
        node_size = []
        colors = mpl.rcParams['axes.prop_cycle'].by_key()['color'] if node_pie_colors is None else node_pie_colors
        if self.node_data is not None:
            for node in T.nodes:
                data = self.node_data.loc[self.node_data['node'] == node].values[0][1:]
                node_size.append(np.sum(data))
                # stat percentage of each node and make the sum of percentage to 1
                pct = np.insert(np.cumsum(data) / np.sum(data), 0, 0)
                # plot scatter-pie
                for i in range(len(pct) - 1):
                    theta = np.linspace(2 * np.pi * pct[i], 2 * np.pi * pct[i + 1], 100)
                    vertices = np.column_stack((np.cos(theta), np.sin(theta)))
                    marker_path = path.Path(np.append(vertices, np.asarray([[0, 0]]), axis=0), closed=False)
                    ax.scatter(*pos[node], np.sum(data) * scale_node_size, colors[i],
                               marker=marker_path, linewidths=0, zorder=2)
            # add node legend
            legend_elements = []
            legend_labels = self.node_data.columns[1:]
            for i, column in enumerate(legend_labels):
                legend_elements.append(
                    patches.Patch(facecolor=colors[i], edgecolor=colors[i], label=column)
                )
            ax.legend(handles=legend_elements,
                      loc='lower left', bbox_to_anchor=(0, 1.02, 1, 0.1),
                      borderaxespad=0, ncols=len(legend_labels), mode='expand', frameon=False)

        # plot the edges with weight
        nx.draw_networkx_edge_labels(
            T, pos, edge_labels={(u, v): d["weight"] * "|" for u, v, d in T.edges(data=True)}
        )
        nx.draw_networkx_edges(T, pos, edge_color=edge_color)
        ax.set_axis_off()
