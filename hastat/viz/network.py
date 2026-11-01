# -*- coding: utf-8 -*-
# @Time    : 2025/7/26 12:58
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : network.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from hastat.viz._export import save_figure
import matplotlib.patches as patches
import networkx as nx
from matplotlib import path
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea, HPacker, TextArea, VPacker
from hastat.log.logger import logger
from hastat.utils.network_io import read_network
from hastat.utils.sample_data import read_sample_csv, sample_frame
from hastat.utils.plot_data import merge_samples


class FancyNetwork(object):
    """
    A class for fancy haplotype network visualization with pie chart nodes
    """
    
    def __init__(self, config):
        """
        Initialize FancyNetwork
        
        :param config: Configuration dictionary containing data and plot parameters
        """
        self.config = config
        self.network = None
        self.sample_hap_data = None
        self.sample_group_data = None
        self.hap_composition = None
        self.node_sizes = None
        self.group_colors = None
        
        # Load and process data
        self._load_data()
        self._get_hap_composition()
        
    def _load_data(self):
        """Load network and sample data"""
        # Load network file
        network_file = self.config['data']['file']
        self.network = read_network(network_file)
        
        logger.info(f"Loaded network with {self.network.number_of_nodes()} nodes and {self.network.number_of_edges()} edges")
        
        # Load sample haplotype data
        sample_hap_file = self.config['data'].get('sample_hap')
        if sample_hap_file:
            self.sample_hap_data = sample_frame(read_sample_csv(sample_hap_file), 'haplotypes').rename(columns={'sample': 'samples'})
            logger.info(f"Loaded sample haplotype data: {len(self.sample_hap_data)} samples")
        
        # Load sample group data
        sample_group_file = self.config['data'].get('sample_group')
        if sample_group_file:
            self.sample_group_data = sample_frame(read_sample_csv(sample_group_file), 'groups').rename(columns={'sample': 'samples'})
            logger.info(f"Loaded sample group data: {len(self.sample_group_data)} samples")
    
    def _get_hap_composition(self):
        """Process haplotype composition data"""
        if self.sample_hap_data is None or self.sample_group_data is None:
            logger.warning("Sample haplotype or group data not provided. Using simple node visualization.")
            self.hap_composition = {}
            self.node_sizes = {node: 1 for node in self.network.nodes()} # default node size is 1
            if self.sample_hap_data is not None:
                self.node_sizes.update(self.sample_hap_data['haplotypes'].value_counts().to_dict())
            return
        
        # Get group index
        group_index = self.config['data'].get('group_index', 1)
        group_columns = [col for col in self.sample_group_data.columns if col != 'samples']
        
        if (not isinstance(group_index, (int, np.integer)) or isinstance(group_index, bool)
                or not 1 <= group_index <= len(group_columns)):
            raise ValueError(f'group_index must be between 1 and {len(group_columns)}')
        
        group_column = group_columns[group_index - 1]
        logger.info(f"Using group column: {group_column}")
        
        # Merge sample data
        merged_data = merge_samples(
            self.sample_hap_data.rename(columns={'samples': 'sample'}),
            self.sample_group_data[['samples', group_column]].rename(columns={'samples': 'sample'})
        )
        
        # Calculate composition for each haplotype
        self.hap_composition = {}
        self.node_sizes = {}
        
        hap_group_counts = merged_data.groupby(['haplotypes', group_column]).size().reset_index(name='count')
        
        for haplotype in self.network.nodes():
            hap_data = hap_group_counts[hap_group_counts['haplotypes'] == haplotype]
            
            if len(hap_data) > 0:
                # Calculate composition percentages
                total_count = hap_data['count'].sum()
                composition = {}
                
                for _, row in hap_data.iterrows():
                    group = row[group_column]
                    count = row['count']
                    composition[group] = count / total_count
                
                self.hap_composition[haplotype] = composition
                self.node_sizes[haplotype] = total_count
            else:
                # No samples for this haplotype
                self.hap_composition[haplotype] = {}
                self.node_sizes[haplotype] = 0
        
        logger.info(f"Processed composition for {len(self.hap_composition)} haplotypes")
        
        # Generate colors for groups
        all_groups = set()
        for comp in self.hap_composition.values():
            all_groups.update(comp.keys())
        
        # Get custom colors from config if provided
        custom_colors = self.config['plot'].get('node_color', None)
        self.group_colors = self._generate_colors(sorted(list(all_groups)), custom_colors)
    
    def _generate_colors(self, groups, custom_colors=None):
        """Generate distinct colors for groups"""
        colors = {}
        n_groups = len(groups)
        
        if n_groups == 0:
            return colors
        cmap = LinearSegmentedColormap.from_list(
            'hastat_groups', ['#C5504B', '#114F8B', '#FCE988', '#90CAEE'], N=100)
        defaults = cmap(np.linspace(0, 1, n_groups))[:, :3]
        
        # Use custom colors if provided
        if custom_colors and len(custom_colors) > 0:
            for i, group in enumerate(groups):
                if i < len(custom_colors):
                    # Use custom color
                    color_str = custom_colors[i]
                    if color_str.startswith('#'):
                        # Convert hex to RGB
                        hex_color = color_str.lstrip('#')
                        rgb = tuple(int(hex_color[i:i+2], 16)/255.0 for i in (0, 2, 4))
                        colors[group] = rgb
                    else:
                        # Try a named color or fall back to the shared palette.
                        try:
                            import matplotlib.colors as mcolors
                            rgb = mcolors.to_rgb(color_str)
                            colors[group] = rgb
                        except:
                            # Keep fallback colors consistent with the other plots.
                            colors[group] = defaults[i]
                else:
                    colors[group] = defaults[i]
        else:
            # Match the palette and sorted-group order used by Pi and other plots.
            for i, group in enumerate(groups):
                colors[group] = defaults[i]
        
        return colors
    
    def _draw_pie_node(self, ax, pos, composition, size, scale_factor=1.0, node_alpha=0.8):
        """
        Draw a pie chart node representing group composition using scatter plots
        
        :param ax: matplotlib axis
        :param pos: node position (x, y)
        :param composition: dictionary of group proportions
        :param size: actual sample count
        :param scale_factor: scaling factor for node sizes
        :param node_alpha: transparency for node colors
        """
        # Calculate node size based on sample count
        node_size = np.sqrt(size) * scale_factor 

        # Mask underlying edges while retaining the requested sector opacity.
        ax.scatter(*pos, s=node_size, c=[ax.get_facecolor()], alpha=1,
                   zorder=1.5, edgecolors='none', linewidth=0)

        # Convert composition to arrays
        groups = sorted(list(composition.keys()))
        proportions = np.array([composition[group] for group in groups])
        colors = [self.group_colors.get(group, 'gray') for group in groups]
        
        if len(groups) == 0:
            # No groups - draw solid circle
            ax.scatter(*pos, s=node_size, c='C3', alpha=node_alpha, 
                zorder=2, edgecolors='white', linewidth=0)
        elif len(groups) == 1:
            # Single group - draw solid circle
            ax.scatter(*pos, s=node_size, c=[colors[0]], alpha=node_alpha, 
                zorder=2, edgecolors='white', linewidth=0)
        else:
            # Multiple groups - draw pie chart using scatter with path markers
            # Calculate cumulative percentages
            pct = np.insert(np.cumsum(proportions), 0, 0)
            
            # Draw each pie slice
            for i in range(len(groups)):
                if proportions[i] > 0:  # Only draw if proportion > 0
                    # Create theta values for this slice
                    theta = np.linspace(2 * np.pi * pct[i], 2 * np.pi * pct[i + 1], 100)
                    
                    # Create vertices for the pie slice
                    vertices = np.column_stack((np.cos(theta), np.sin(theta)))
                    # Add center point to close the slice
                    vertices = np.append(vertices, np.array([[0, 0]]), axis=0)
                    
                    # Create path for this slice
                    marker_path = path.Path(vertices, closed=False)
                    
                    # Draw the slice
                    ax.scatter(*pos, s=node_size, c=[colors[i]], 
                             marker=marker_path, alpha=node_alpha, 
                             zorder=2, edgecolors='none', linewidth=0)
    
    def _get_layout_positions(self, layout='spring', seed=None, k=None):
        """Get node positions using specified layout algorithm"""
        layout_funcs = {
            'spring': nx.spring_layout,
            'circular': nx.circular_layout,
            'kamada_kawai': nx.kamada_kawai_layout,
            'fruchterman_reingold': nx.fruchterman_reingold_layout,
            'shell': nx.shell_layout,
            'spectral': nx.spectral_layout,
            'random': nx.random_layout
        }
        
        if layout not in layout_funcs:
            logger.warning(f"Unknown layout '{layout}', using 'spring' layout")
            layout = 'spring'
        
        logger.info(f"Using {layout} layout" + (f" with seed {seed}" if seed is not None else ""))
        
        try:
            if layout == 'kamada_kawai' and self.network.number_of_nodes() > 1:
                pos = layout_funcs[layout](self.network, weight='weight')
            elif layout == 'spring':
                pos = layout_funcs[layout](self.network, seed=seed, k=k)
            elif layout in ['fruchterman_reingold', 'random']:
                # These layouts support seed parameter
                pos = layout_funcs[layout](self.network, seed=seed)
            else:
                pos = layout_funcs[layout](self.network)
        except:
            logger.warning(f"Failed to use {layout} layout, falling back to spring layout")
            pos = nx.spring_layout(self.network, seed=seed, k=k)
        
        return pos
    
    def _draw_curved_edge(self, ax, pos1, pos2, edge_color='gray', edge_width=1.0, edge_alpha=0.6, curvature=0.2):
        """
        Draw a curved edge between two nodes
        
        :param ax: matplotlib axis
        :param pos1: position of first node (x1, y1)
        :param pos2: position of second node (x2, y2)
        :param edge_color: edge color
        :param edge_width: edge width
        :param edge_alpha: edge transparency
        :param curvature: curvature factor (0 = straight line, higher = more curved)
        """
        x1, y1 = pos1
        x2, y2 = pos2
        
        # Calculate control point for Bezier curve
        mid_x = (x1 + x2) / 2
        mid_y = (y1 + y2) / 2
        
        # Calculate perpendicular vector for curvature
        dx = x2 - x1
        dy = y2 - y1
        length = np.sqrt(dx**2 + dy**2)
        
        if length > 0:
            # Perpendicular vector (rotated 90 degrees)
            perp_x = -dy / length
            perp_y = dx / length
            
            # Control point offset by curvature
            control_x = mid_x + perp_x * curvature * length
            control_y = mid_y + perp_y * curvature * length
            
            # Create Bezier curve points
            t = np.linspace(0, 1, 50)
            bezier_x = (1-t)**2 * x1 + 2*(1-t)*t * control_x + t**2 * x2
            bezier_y = (1-t)**2 * y1 + 2*(1-t)*t * control_y + t**2 * y2
            
            ax.plot(bezier_x, bezier_y, color=edge_color, 
                   linewidth=edge_width, alpha=edge_alpha, zorder=1)
        else:
            # Fallback to straight line if nodes are at same position
            ax.plot([x1, x2], [y1, y2], color=edge_color, 
                   linewidth=edge_width, alpha=edge_alpha, zorder=1)
    
    def plot(self):
        """Create the network visualization"""
        import matplotlib
        matplotlib.use('Agg')  # Set backend before importing pyplot
        
        plot_config = self.config['plot']
        
        # Set up figure
        width = plot_config.get('width', 10)
        height = plot_config.get('height', 8) 
        fig, ax = plt.subplots(figsize=(width, height))
        
        # Get layout
        layout = plot_config.get('layout', 'spring')
        seed = plot_config.get('seed', None)
        k = plot_config.get('k', None)
        pos = self._get_layout_positions(layout, seed, k)
        
        # Draw edges
        edge_color = plot_config.get('edge_color', 'gray')
        edge_width = plot_config.get('edge_width', 1.0)
        edge_alpha = plot_config.get('edge_alpha', 0.6)
        edge_style = plot_config.get('edge_style', 'line')
        
        if edge_style == 'curve':
            # Use custom curved edge drawing
            for u, v, data in self.network.edges(data=True):
                pos1 = pos[u]
                pos2 = pos[v]
                self._draw_curved_edge(ax, pos1, pos2, edge_color, edge_width, edge_alpha)
        else:
            # Use NetworkX built-in edge drawing for straight lines
            nx.draw_networkx_edges(self.network, pos, edge_color=edge_color, 
                                 width=edge_width, alpha=edge_alpha, ax=ax)
        
        # plot node
        scale_factor = plot_config.get('node_scale_factor', 1.0)
        
        # Get node alpha
        node_alpha = plot_config.get('node_alpha', 0.8)
        
        # Draw nodes and collect scatter objects for legend
        for node in self.network.nodes():
            node_pos = pos[node]
            composition = self.hap_composition.get(node, {})
            size = self.node_sizes.get(node, 1)
            
            self._draw_pie_node(ax, node_pos, composition, size, 
                scale_factor=scale_factor, node_alpha=node_alpha)
        
        # Draw node labels if requested
        show_node_label = plot_config.get('show_node_label', False)
        node_label_size = plot_config.get('node_label_size', 10)
        node_label_color = plot_config.get('node_label_color', 'black')
        
        if show_node_label:
            nx.draw_networkx_labels(self.network, pos, font_size=node_label_size, 
                                  font_weight='bold', font_color=node_label_color, ax=ax)
        
        # Draw edge labels if requested
        show_edge_label = plot_config.get('show_edge_label', True)
        edge_label_symbol = plot_config.get('edge_label_symbol', None)
        edge_label_size = plot_config.get('edge_label_size', 8)
        edge_label_color = plot_config.get('edge_label_color', 'black')
        
        if show_edge_label:
            if edge_label_symbol:
                # Create edge labels with symbol multiplication (weight * symbol)
                edge_labels = {(u, v): edge_label_symbol * int(d.get('weight', 1)) 
                             for u, v, d in self.network.edges(data=True)}
            else:
                # Show plain weight numbers
                edge_labels = {(u, v): str(d.get('weight', 1)) 
                             for u, v, d in self.network.edges(data=True)}
            
            nx.draw_networkx_edge_labels(self.network, pos, edge_labels=edge_labels,
                                       font_size=edge_label_size, font_color=edge_label_color, ax=ax)
        
        # Create legends
        x, y, node_size = [], [], []
        for node in self.network.nodes():
            node_pos = pos[node]
            x.append(node_pos[0])
            y.append(node_pos[1])
            node_size.append(np.sqrt(self.node_sizes.get(node, 1)) * scale_factor)
        scatter = ax.scatter(x, y, s=node_size, c='none', zorder=2, edgecolors='white', linewidth=0.5)
        self._create_legends(ax, scatter)
        
        # Set axis properties
        ax.axis('off')
        
        # Set title
        title = plot_config.get('title', 'Haplotype Network')
        ax.set_title(title, fontsize=14, weight='bold', pad=20)
        
        # Save figure
        save_fig = plot_config.get('save_fig')
        if save_fig:
            save_figure(fig, save_fig, dpi=300, bbox_inches='tight', facecolor='white')
            logger.info(f"Network plot saved to: {save_fig}")
        
        # Don't call plt.show() in headless environment
        plt.close(fig)
    
    def _create_legends(self, ax, scatter):
        """Pack readable legends outside the network's plotting area."""
        sections = []
        font_size = 10
        if self.group_colors:
            rows = [TextArea('Population', textprops=dict(size=font_size, weight='bold'))]
            for group, color in self.group_colors.items():
                swatch = DrawingArea(12, 10, 0, 0)
                swatch.add_artist(patches.Rectangle((0, 0), 12, 10, facecolor=color,
                    linewidth=0, alpha=self.config['plot'].get('node_alpha', .8)))
                rows.append(HPacker(children=[swatch, TextArea(str(group), textprops=dict(size=font_size))],
                                    pad=0, sep=6, align='center'))
            sections.append(VPacker(children=rows, pad=0, sep=5, align='left'))
        sizes = self._create_size_legend(ax, scatter)
        if sizes is not None:
            sections.append(sizes)
        if sections:
            self.legend_box = AnchoredOffsetbox(loc='upper left', pad=0, borderpad=0,
                child=VPacker(children=sections, pad=0, sep=16, align='left'),
                bbox_to_anchor=(1.03, 1), bbox_transform=ax.transAxes, frameon=False)
            self.legend_box.set_clip_on(False)
            ax.add_artist(self.legend_box)

    def _size_legend_counts(self):
        """Use actual minimum, middle and maximum nonzero sample counts."""
        if self.sample_hap_data is None:
            return []
        counts = sorted({int(self.node_sizes.get(node, 0)) for node in self.network
                         if self.node_sizes.get(node, 0) > 0})
        if len(counts) <= 3:
            return counts
        return [counts[0], counts[len(counts)//2], counts[-1]]
    
    def _create_size_legend(self, ax, scatter):
        """Keep legend circles at exactly the same point sizes as network nodes."""
        scale_factor = self.config['plot'].get('node_scale_factor', 1.0)
        label_spacing = self.config['plot'].get('legend_label_spacing', 2)
        handle_text_pad = self.config['plot'].get('legend_handle_text_pad', 1)
        counts = self._size_legend_counts()
        if not counts:
            return None
        font_size = 10
        diameters = np.sqrt(np.sqrt(counts) * scale_factor)
        width = max(diameters)
        rows = []
        self.size_legend_markers = []
        self.size_legend_labels = []
        for count, diameter in zip(counts, diameters):
            handle = DrawingArea(width + 2, diameter + 2, 0, 0)
            circle = patches.Circle(((width + 2)/2, (diameter + 2)/2), diameter/2,
                                    facecolor=ax.get_facecolor(), edgecolor='0.4', linewidth=.7)
            handle.add_artist(circle)
            label = TextArea(f'{count:,}', textprops=dict(size=font_size))
            rows.append(HPacker(children=[handle, label], pad=0,
                                sep=max(0, handle_text_pad)*font_size, align='center'))
            self.size_legend_markers.append(circle)
            self.size_legend_labels.append(label)
        return VPacker(children=[TextArea('Samples (n)', textprops=dict(size=font_size, weight='bold')),
            VPacker(children=rows, pad=0, sep=max(0, label_spacing)*font_size, align='left')],
            pad=0, sep=8, align='left')
