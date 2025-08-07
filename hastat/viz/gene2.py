# -*- coding: utf-8 -*-
# @Time    : 2025/07/25 12:58
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : gene2.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import tomli
from hastat.log.logger import logger
from hastat.dataio.gff import read as read_gff
from hastat.utils.gene import GeneFeature


class FancyGene(object):
    """
    A class for fancy gene visualization using matplotlib based on TOML configuration files.
    
    This class extracts gene structures (CDS, exons, UTRs) from GFF files and visualizes them
    according to layout requirements specified in TOML configuration files. It can also
    optionally visualize annotation information and haplotype data.
    """
    
    def __init__(self, config: dict = None):
        """
        Initialize FancyGene object.
        
        :param config: Configuration dictionary with 'data' and 'plot' sections
        """
        if config is None:
            raise ValueError("config is None")
        self.config = config
        
        # Data attributes - will be set by add_data()
        self.gff_file = None
        self.genes = None
        self.toml_file = None
        self.toml_config = None
        self.gff_db = None
        self.gene_feature = None
        
        # Plot attributes - from config
        self.upstream = config['plot'].get('upstream', 0)
        self.downstream = config['plot'].get('downstream', 0)
        self.width = config['plot'].get('width', 12.0)
        self.height = config['plot'].get('height', 8.0)
        self.output = config['plot'].get('save_fig', None)
        
        # Gene data storage
        self.gene_data = {}
        self.gene_regions = {}
        
        # Initialize figure and axes
        self.fig = None
        self.ax = None
        
    def add_data(self, data: dict):
        """
        Add data from the config data section.
        
        :param data: Data dictionary containing gff_file, genes, toml_file
        """
        self.gff_file = data['gff_file']
        self.genes = data['genes'] if isinstance(data['genes'], list) else [data['genes']]
        self.toml_file = data['toml_file']
        
        # Load TOML configuration
        self.toml_config = self._load_toml_config()
        
        # Load GFF database
        self.gff_db = read_gff(self.gff_file)
        self.gene_feature = GeneFeature(self.gff_db)
        
    def _load_toml_config(self):
        """Load TOML configuration file."""
        try:
            with open(self.toml_file, 'rb') as f:
                config = tomli.load(f)
            logger.info(f"Successfully loaded TOML configuration: {self.toml_file}")
            return config
        except Exception as e:
            logger.error(f"Failed to load TOML configuration: {e}")
            raise e
    
    def _get_gene_config(self, gene_id):
        """Get configuration for a specific gene, falling back to default if not found."""
        if gene_id in self.toml_config:
            # Merge with default config
            gene_config = self.toml_config['default'].copy()
            gene_config.update(self.toml_config[gene_id])
            return gene_config
        else:
            logger.warning(f"Gene {gene_id} not found in config, using default values")
            return self.toml_config['default'].copy()
    
    def _extract_gene_features(self, gene_id):
        """Extract gene features (CDS, exons, UTRs) from GFF database."""
        try:
            # Get gene locus with upstream/downstream regions
            chrom, start, end = self.gene_feature.get_locus(
                gene_id, upstream=self.upstream, downstream=self.downstream
            )
            
            # Get the gene object
            gene_obj = self.gff_db[gene_id]
            
            # Extract features
            features = {
                'gene_id': gene_id,
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': gene_obj.strand,
                'gene_start': gene_obj.start,
                'gene_end': gene_obj.end,
                'CDS': [],
                'exon': [],
                'five_prime_UTR': [],
                'three_prime_UTR': [],
                'UTR': []
            }
            
            # Get all features for this gene
            for feature in self.gff_db.children(gene_obj, order_by='start'):
                feature_type = feature.featuretype
                if feature_type in features:
                    features[feature_type].append([feature.start, feature.end])
            
            # Combine UTRs
            features['UTR'] = features['five_prime_UTR'] + features['three_prime_UTR']
            
            logger.info(f"Extracted features for gene {gene_id}: {chrom}:{start}-{end}")
            return features
            
        except Exception as e:
            logger.error(f"Failed to extract features for gene {gene_id}: {e}")
            raise e
    
    def _load_annotation_data(self, ann_file):
        """Load annotation data from CSV file."""
        if not ann_file or ann_file == "":
            return None
        try:
            ann_data = pd.read_csv(ann_file)
            logger.info(f"Loaded annotation data from {ann_file}")
            # remove nan and empty string
            if ann_data is not None:
                ann_data = ann_data.dropna()
                ann_data = ann_data[ann_data['ann'] != '']

            logger.info(f"detected {len(ann_data)} variant annotations")
            return ann_data
        except Exception as e:
            logger.warning(f"Failed to load annotation file {ann_file}: {e}")
            return None
    
    def _load_haplotype_data(self, hap_file, hap_cols=None):
        """Load haplotype data from CSV file."""
        if not hap_file or hap_file == "":
            return None
        try:
            hap_data = pd.read_csv(hap_file)
            logger.info(f"Loaded haplotype data from {hap_file}")
            if isinstance(hap_cols, list) and len(hap_cols) > 0:
                hap_data = hap_data.loc[:, ['chrom', 'pos', 'ref', 'alt'] + hap_cols]
            logger.info(f"detected {len(hap_data.columns) - 4} haplotypes with {len(hap_data)} variant sites")
            return hap_data
        except Exception as e:
            logger.warning(f"Failed to load haplotype file {hap_file}: {e}")
            return None
    
    def _coordinate_transform(self, pos, bbox, gene_start, gene_end):
        """Transform genomic coordinates to plot coordinates based on bbox."""
        x_start, y_start, width, height = bbox
        # Transform genomic position to relative position within gene
        relative_pos = (pos - gene_start) / (gene_end - gene_start)
        # Transform to plot coordinates
        plot_x = x_start + relative_pos * width
        return plot_x
    
    def _draw_gene_structure(self, gene_data, gene_config, bbox):
        """Draw gene structure features (CDS, exons, UTRs)."""
        x_start, y_start, bbox_width, bbox_height = bbox
        xy = gene_config.get('xy', [0, 0.5])
        width = gene_config.get('width', 1.0)
        height = gene_config.get('height', 0.2)
        color = gene_config.get('color', 'C1')
        feature_type = gene_config.get('feature', 'CDS')
        
        # Calculate absolute coordinates
        abs_x = x_start + xy[0] * bbox_width
        abs_y = y_start + xy[1] * bbox_height
        abs_width = width * bbox_width
        abs_height = height * bbox_height
        
        gene_start = gene_data['gene_start']
        gene_end = gene_data['gene_end']
        
        # Draw selected features
        features_to_draw = gene_data.get(feature_type, [])
        if len(features_to_draw) == 0:
            logger.warning(f"No {feature_type} features found! Please check the GFF file and toml file.")
        for feature_start, feature_end in features_to_draw:
            # Calculate relative positions
            rel_start = (feature_start - gene_start) / (gene_end - gene_start)
            rel_end = (feature_end - gene_start) / (gene_end - gene_start)
            
            # Draw feature rectangle
            feature_x = abs_x + rel_start * abs_width
            feature_width = (rel_end - rel_start) * abs_width
            
            feature_rect = patches.Rectangle(
                (feature_x, abs_y - abs_height/2),
                feature_width, abs_height,
                facecolor=color, edgecolor=color, 
                alpha=1, zorder=2
            )
            self.ax.add_patch(feature_rect)
        
        # Draw strand arrow if gene has features
        if len(features_to_draw) > 0:
            arrow_y = abs_y
            arrow_extension = abs_width * 0.05  # arrow extension length: 5% of the gene width
            
            if gene_data['strand'] == '+':
                # + strand: arrow from the start of the gene to the end of the gene
                arrow_start = abs_x - arrow_extension
                arrow_end = abs_x + abs_width + arrow_extension
                arrow = patches.FancyArrowPatch(
                    (arrow_start, arrow_y), (arrow_end, arrow_y),
                    arrowstyle='->', mutation_scale=15, color='black', 
                    alpha=0.6, zorder=1
                )
            else:
                # - strand: arrow from the end of the gene to the start of the gene
                arrow_start = abs_x + abs_width + arrow_extension
                arrow_end = abs_x - arrow_extension
                arrow = patches.FancyArrowPatch(
                    (arrow_start, arrow_y), (arrow_end, arrow_y),
                    arrowstyle='->', mutation_scale=15, color='black', 
                    alpha=0.6, zorder=1
                )
            self.ax.add_patch(arrow)
        
        return abs_x, abs_y, abs_width, abs_height
    
    def _draw_annotations(self, gene_data, gene_config, bbox, gene_coords):
        """Draw variant annotations."""
        anno_colors = ['#C5504B', '#114F8B', '#FCE988', '#90CAEE']  # Distinct colors
        ann_file = gene_config.get('ann_file', '')
        ann_y = gene_config.get('ann_y', 0.8)
        ann_text_size = gene_config.get('ann_text_size', 8)
        
        ann_data = self._load_annotation_data(ann_file)

        if ann_data is None or ann_data.empty:
            return

        x_start, y_start, bbox_width, bbox_height = bbox
        abs_x, abs_y, abs_width, abs_height = gene_coords
        
        # Filter annotations for this gene region
        gene_start = gene_data['gene_start']
        gene_end = gene_data['gene_end']
        chrom = gene_data['chrom']
        
        # Assuming annotation data has columns: chrom, pos, ann
        if 'chrom' in ann_data.columns:
            gene_ann = ann_data[
                (ann_data['chrom'] == chrom) & 
                (ann_data['pos'] >= gene_start) & 
                (ann_data['pos'] <= gene_end)
            ]
        else:
            gene_ann = ann_data[
                (ann_data['pos'] >= gene_start) & 
                (ann_data['pos'] <= gene_end)
            ]
        
        if gene_ann.empty:
            return
        
        # Calculate original positions and prepare annotation data
        annotations = []
        ann_abs_y = y_start + ann_y * bbox_height
        
        for i, (_, row) in enumerate(gene_ann.iterrows()):
            pos = row['pos']
            annotation = row['ann']
            color = anno_colors[i % len(anno_colors)]
            
            # Calculate x position on gene
            rel_pos = (pos - gene_start) / (gene_end - gene_start)
            original_x = abs_x + rel_pos * abs_width
            
            annotations.append({
                'pos': pos,
                'original_x': original_x,
                'text': annotation,
                'color': color,
                'width': len(annotation) * 0.008 * abs_width  # Estimate text width
            })
        
        # Sort annotations by original x position
        annotations.sort(key=lambda x: x['original_x'])
        
        # Resolve text overlaps using a simple greedy algorithm
        adjusted_annotations = self._resolve_annotation_overlaps(annotations, abs_x, abs_width)
        
        # Draw annotation markers with adjusted positions
        for ann in adjusted_annotations:
            # Text at adjusted position at ann_y height
            text_x = ann['adjusted_x']
            text_y = ann_abs_y
            
            # Arrow target: original position on gene
            arrow_target_x = ann['original_x']
            arrow_target_y = abs_y + abs_height/2  # Top of gene features
            
            # Draw annotation marker with arrow pointing from text to gene position
            self.ax.annotate(
                ann['text'], 
                xy=(arrow_target_x, arrow_target_y),  # Arrow points to this position on gene
                xytext=(text_x, text_y),  # Text is at this position
                ha='center', va='bottom', rotation='vertical', fontsize=ann_text_size,
                color=ann['color'], weight='bold',
                bbox=dict(boxstyle='round,pad=0.3', 
                        facecolor='white', 
                        edgecolor=ann['color'], 
                        alpha=0.8),
                arrowprops=dict(
                        arrowstyle='->', 
                        color=ann['color'], 
                        alpha=0.7,
                        connectionstyle="arc3,rad=0")
            )
    
    def _resolve_annotation_overlaps(self, annotations, gene_x_start, gene_width):
        """
        Resolve text overlaps in annotation layout using a greedy algorithm.
        
        :param annotations: List of annotation dictionaries
        :param gene_x_start: Start x position of gene
        :param gene_width: Width of gene region
        :return: List of annotations with adjusted x positions
        """
        if not annotations:
            return []
        
        # Start with original positions
        for ann in annotations:
            ann['adjusted_x'] = ann['original_x']
        
        # Iteratively resolve overlaps
        changed = True
        max_iterations = 10
        iteration = 0
        
        while changed and iteration < max_iterations:
            changed = False
            iteration += 1
            
            for i in range(len(annotations) - 1):
                curr = annotations[i]
                next_ann = annotations[i + 1]
                
                # Check for overlap (considering text width)
                overlap = (curr['adjusted_x'] + curr['width']/2) - (next_ann['adjusted_x'] - next_ann['width']/2)
                
                if overlap > 0:
                    # Calculate adjustment needed
                    adjustment = overlap / 2 + 0.01  # Small buffer
                    
                    # Move current annotation left and next annotation right
                    new_curr_x = curr['adjusted_x'] - adjustment
                    new_next_x = next_ann['adjusted_x'] + adjustment
                    
                    # Ensure annotations stay within reasonable bounds
                    min_x = gene_x_start - gene_width * 0.3
                    max_x = gene_x_start + gene_width * 1.3
                    
                    if new_curr_x >= min_x:
                        curr['adjusted_x'] = new_curr_x
                        changed = True
                    
                    if new_next_x <= max_x:
                        next_ann['adjusted_x'] = new_next_x
                        changed = True
        
        return annotations
    
    def _draw_haplotypes(self, gene_data, gene_config, bbox):
        """Draw haplotype information."""
        hap_file = gene_config.get('hap_file', '')
        hap_y = gene_config.get('hap_y', 0)
        hap_height = gene_config.get('hap_height', 0.3)
        hap_label_show = gene_config.get('hap_label_show', True)
        hap_label_size = gene_config.get('hap_label_size', 8)
        hap_cols = gene_config.get('hap_cols', [])
        hap_data = self._load_haplotype_data(hap_file, hap_cols)
        if hap_data is None:
            return
        
        x_start, y_start, bbox_width, bbox_height = bbox
        
        # Filter haplotype data for this gene region
        chrom_start = gene_data['start']
        chrom_end = gene_data['end']
        chrom = gene_data['chrom']
        
        if 'chrom' in hap_data.columns:
            gene_hap = hap_data[
                (hap_data['chrom'] == chrom) & 
                (hap_data['pos'] >= chrom_start) & 
                (hap_data['pos'] <= chrom_end)
            ].copy()
        else:
            gene_hap = hap_data[
                (hap_data['pos'] >= chrom_start) & 
                (hap_data['pos'] <= chrom_end)
            ].copy()
        
        if gene_hap.empty:
            return
        
        # Get haplotype columns (excluding chrom, pos, ref, alt)
        hap_cols = [col for col in gene_hap.columns 
                   if col not in ['chrom', 'pos', 'ref', 'alt']]
        
        if not hap_cols:
            return

        # Calculate coordinates for equal-spaced grid
        hap_abs_x = x_start
        hap_abs_y = y_start + hap_y * bbox_height
        hap_abs_width = bbox_width
        hap_abs_height = hap_height * bbox_height
        
        num_haps = len(hap_cols)  # Number of rows (haplotypes)
        num_sites = len(gene_hap)  # Number of columns (variant sites)
        
        # Calculate cell dimensions
        cell_width = hap_abs_width / num_sites
        cell_height = hap_abs_height / num_haps
        
        # Draw haplotype matrix with equal spacing
        colors = {0: 'white', 1: '#90CAEE', 2: '#114F8B', -1: 'lightgray'}
        
        # Sort gene_hap by position to ensure consistent column order
        gene_hap_sorted = gene_hap.sort_values('pos')
        
        for hap_idx, hap_col in enumerate(hap_cols):
            for site_idx, (_, row) in enumerate(gene_hap_sorted.iterrows()):
                genotype = row[hap_col]
                
                # Calculate cell position (equal-spaced grid)
                cell_x = hap_abs_x + site_idx * cell_width
                cell_y = hap_abs_y + hap_idx * cell_height
                
                # Draw genotype rectangle
                color = colors.get(genotype, 'gray')
                rect = patches.Rectangle(
                    (cell_x, cell_y),
                    cell_width, cell_height,
                    facecolor=color, edgecolor='black', linewidth=0.5
                )
                self.ax.add_patch(rect)
        
        # Add haplotype labels on the left
        if hap_label_show:
            for hap_idx, hap_col in enumerate(hap_cols):
                label_y = hap_abs_y + hap_idx * cell_height + cell_height/2
                self.ax.text(
                    hap_abs_x - 0.02, label_y,
                    hap_col, ha='right', va='center', fontsize=hap_label_size
                )
        
        # Draw connecting lines from haplotype columns to gene variant positions
        self._draw_haplotype_connections(gene_data, gene_hap_sorted, bbox, 
                                       hap_abs_x, hap_abs_y, hap_abs_height, 
                                       cell_width)
    
    def _draw_haplotype_connections(self, gene_data, gene_hap_sorted, bbox, 
                                  hap_abs_x, hap_abs_y, hap_abs_height, cell_width):
        """Draw connecting lines from haplotype columns to gene variant positions."""
        x_start, y_start, bbox_width, bbox_height = bbox
        gene_start = gene_data['gene_start']
        gene_end = gene_data['gene_end']
        
        # Calculate gene region coordinates (same as in _draw_gene_structure)
        gene_config = self._get_gene_config(gene_data['gene_id'])
        xy = gene_config.get('xy', [0, 0.5])
        width = gene_config.get('width', 1.0)
        height = gene_config.get('height', 0.2)
        
        # Calculate absolute coordinates for gene
        gene_abs_x = x_start + xy[0] * bbox_width
        gene_abs_y = y_start + xy[1] * bbox_height
        gene_abs_width = width * bbox_width
        gene_abs_height = height * bbox_height
        
        # Draw connection lines for each variant site
        for site_idx, (_, row) in enumerate(gene_hap_sorted.iterrows()):
            pos = row['pos']
            
            # Position in haplotype matrix (top of the column)
            hap_column_x = hap_abs_x + site_idx * cell_width + cell_width/2
            hap_column_y_top = hap_abs_y + hap_abs_height  # Top of haplotype matrix
            
            # Position on gene (variant location)
            rel_pos = (pos - gene_start) / (gene_end - gene_start)
            gene_variant_x = gene_abs_x + rel_pos * gene_abs_width
            gene_variant_y = gene_abs_y - gene_abs_height  # Bottom of gene
            
            # Draw connecting line: diagonal line first, then vertical line to gene bottom
            # Connection height is half the distance between haplotype top and gene bottom
            connection_height = abs(gene_variant_y - hap_column_y_top) * 0.5
            intermediate_y = gene_variant_y - connection_height

            # Diagonal line from haplotype column to intermediate point
            self.ax.plot([hap_column_x, gene_variant_x], 
                        [hap_column_y_top, intermediate_y],
                        color='gray', alpha=0.6, linewidth=1, linestyle='-')
            
            # Vertical line up to gene variant position (bottom of gene)
            self.ax.plot([gene_variant_x, gene_variant_x], 
                        [intermediate_y, gene_variant_y],
                        color='gray', alpha=0.6, linewidth=1, linestyle='-')

    def _load_pi_data(self, pi_file):
        """Load pi data from CSV file."""
        if not pi_file or pi_file == "":
            return None
        try:
            pi_data = pd.read_csv(pi_file)
            logger.info(f"Loaded pi data from {pi_file}")
            logger.info(f"Detected {len(pi_data)} pi data points")
            return pi_data
        except Exception as e:
            logger.warning(f"Failed to load pi file {pi_file}: {e}")
            return None
    
    def _load_fst_data(self, fst_file):
        """Load fst data from CSV file."""
        if not fst_file or fst_file == "":
            return None
        try:
            fst_data = pd.read_csv(fst_file)
            logger.info(f"Loaded fst data from {fst_file}")
            logger.info(f"Detected {len(fst_data)} fst data points")
            return fst_data
        except Exception as e:
            logger.warning(f"Failed to load fst file {fst_file}: {e}")
            return None
    
    def _draw_pi_plot(self, gene_data, gene_config, bbox):
        """Draw pi data plot using inset axes."""
        pi_file = gene_config.get('pi_file', '')
        pi_y = gene_config.get('pi_y', 0.55)
        pi_height = gene_config.get('pi_height', 0.2)
        pi_style = gene_config.get('pi_style', 'fill')
        pi_cmap = LinearSegmentedColormap.from_list('pi_cmap', ['#C5504B', '#114F8B', '#FCE988', '#90CAEE'], N=100)
        
        pi_data = self._load_pi_data(pi_file)
        if pi_data is None or pi_data.empty:
            return None
        
        x_start, y_start, bbox_width, bbox_height = bbox
        
        # Filter pi data for this gene region
        chrom_start = gene_data['start']
        chrom_end = gene_data['end']
        chrom = gene_data['chrom']
        
        if 'chrom' in pi_data.columns:
            gene_pi = pi_data[
                (pi_data['chrom'] == chrom) & 
                (pi_data['start'] >= chrom_start) & 
                (pi_data['end'] <= chrom_end)
            ].copy()
        else:
            gene_pi = pi_data[
                (pi_data['start'] >= chrom_start) & 
                (pi_data['end'] <= chrom_end)
            ].copy()
        
        if gene_pi.empty:
            return None
        
        # Calculate genomic positions (center of each window)
        gene_pi['pos_bp'] = (gene_pi['start'] + gene_pi['end']) / 2
        
        # Convert pi values to scientific notation (multiply by 1000 for 10^-3 scale)
        gene_pi['pi_scaled'] = gene_pi['pi'] * 1000
        
        # Create inset axes for pi plot
        pi_bounds = [x_start, y_start + pi_y * bbox_height, bbox_width, pi_height * bbox_height]
        pi_ax = self.ax.inset_axes(pi_bounds, transform=self.ax.transData)
        
        # Draw plots based on style and groups
        legend_handles = []
        
        if 'groups' in gene_pi.columns:
            groups = sorted(gene_pi['groups'].unique())
            colors = pi_cmap(np.linspace(0, 1, len(groups)))
            
            for i, group in enumerate(groups):
                group_data = gene_pi[gene_pi['groups'] == group].sort_values('pos_bp')
                if group_data.empty:
                    continue
                
                color = colors[i]
                
                if pi_style == 'fill':
                    # Fill between style
                    pi_ax.fill_between(group_data['pos_bp'], 0, group_data['pi_scaled'], 
                                     color=color, alpha=0.6, label=group)
                    line = pi_ax.plot(group_data['pos_bp'], group_data['pi_scaled'], 
                                    color=color, linewidth=1, alpha=0.8, label=group)[0]
                else:
                    # Line only style
                    line = pi_ax.plot(group_data['pos_bp'], group_data['pi_scaled'], 
                                    color=color, linewidth=1, label=group)[0]
                
                legend_handles.append(line)
        else:
            # Single group
            sorted_data = gene_pi.sort_values('pos_bp')
            color = 'C0'
            
            if pi_style == 'fill':
                pi_ax.fill_between(sorted_data['pos_bp'], 0, sorted_data['pi_scaled'], 
                                 color=color, alpha=0.6, label='π')
                line = pi_ax.plot(sorted_data['pos_bp'], sorted_data['pi_scaled'], 
                                color=color, linewidth=1, alpha=0.8, label='π')[0]
            else:
                line = pi_ax.plot(sorted_data['pos_bp'], sorted_data['pi_scaled'], 
                                color=color, linewidth=1, label='π')[0]
            
            legend_handles = [line]
        
        # Configure pi axes
        x_min, x_max = np.min(gene_pi['pos_bp']), np.max(gene_pi['pos_bp'])
        pi_ax.set_xlim(x_min, x_max)
        pi_ax.set_ylabel(r'π ($10^{-3}$)', fontsize=8)
        
        # Set custom x-axis ticks only at min and max values
        pi_ax.set_xticks([x_min, x_max])
        pi_ax.set_xticklabels([f'{int(x_min):,}', f'{int(x_max):,}'])
        pi_ax.tick_params(axis='both', labelsize=6)
        pi_ax.spines[["top", "right"]].set_visible(False)
        
        # Add legend
        if legend_handles:
            pi_ax.legend(legend_handles, [h.get_label() for h in legend_handles],
                        bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                        mode='expand', borderaxespad=0., ncol=3, fontsize=6, frameon=False)
        
        return legend_handles
    
    def _draw_fst_plot(self, gene_data, gene_config, bbox):
        """Draw fst data plot using inset axes."""
        fst_file = gene_config.get('fst_file', '')
        fst_y = gene_config.get('fst_y', 0.8)
        fst_height = gene_config.get('fst_height', 0.2)
        fst_style = gene_config.get('fst_style', 'fill')
        fst_cmap = LinearSegmentedColormap.from_list('fst_cmap', ['#C5504B', '#114F8B', '#FCE988', '#90CAEE'], N=100)

        fst_data = self._load_fst_data(fst_file)
        if fst_data is None or fst_data.empty:
            return None
        
        x_start, y_start, bbox_width, bbox_height = bbox
        
        # Filter fst data for this gene region
        chrom_start = gene_data['start']
        chrom_end = gene_data['end']
        chrom = gene_data['chrom']
        
        if 'chrom' in fst_data.columns:
            gene_fst = fst_data[
                (fst_data['chrom'] == chrom) & 
                (fst_data['start'] >= chrom_start) & 
                (fst_data['end'] <= chrom_end)
            ].copy()
        else:
            gene_fst = fst_data[
                (fst_data['start'] >= chrom_start) & 
                (fst_data['end'] <= chrom_end)
            ].copy()
        
        if gene_fst.empty:
            return None
        
        # Calculate genomic positions (center of each window)
        gene_fst['pos_bp'] = (gene_fst['start'] + gene_fst['end']) / 2
        
        # Create inset axes for fst plot
        fst_bounds = [x_start, y_start + fst_y * bbox_height, bbox_width, fst_height * bbox_height]
        fst_ax = self.ax.inset_axes(fst_bounds, transform=self.ax.transData)
        
        # Draw plots based on style and groups
        legend_handles = []
        
        if 'groups' in gene_fst.columns:
            groups = sorted(gene_fst['groups'].unique())
            colors = fst_cmap(np.linspace(0, 1, len(groups)))
            
            for i, group in enumerate(groups):
                group_data = gene_fst[gene_fst['groups'] == group].sort_values('pos_bp')
                if group_data.empty:
                    continue
                
                color = colors[i]
                
                if fst_style == 'fill':
                    # Fill between style
                    fst_ax.fill_between(group_data['pos_bp'], 0, group_data['fst'], 
                                      color=color, alpha=0.6, label=group)
                    line = fst_ax.plot(group_data['pos_bp'], group_data['fst'], 
                                     color=color, linewidth=1, alpha=0.8, label=group)[0]
                else:
                    # Line only style
                    line = fst_ax.plot(group_data['pos_bp'], group_data['fst'], 
                                     color=color, linewidth=1, label=group)[0]
                
                legend_handles.append(line)
        else:
            # Single group
            sorted_data = gene_fst.sort_values('pos_bp')
            color = 'C1'
            
            if fst_style == 'fill':
                fst_ax.fill_between(sorted_data['pos_bp'], 0, sorted_data['fst'], 
                                  color=color, alpha=0.6, label='F$_{st}$')
                line = fst_ax.plot(sorted_data['pos_bp'], sorted_data['fst'], 
                                 color=color, linewidth=1.5, alpha=0.8, label='F$_{st}$')[0]
            else:
                line = fst_ax.plot(sorted_data['pos_bp'], sorted_data['fst'], 
                                 color=color, linewidth=2, label='F$_{st}$')[0]
            
            legend_handles = [line]
        
        # Configure fst axes
        x_min, x_max = np.min(gene_fst['pos_bp']), np.max(gene_fst['pos_bp'])
        fst_ax.set_xlim(x_min, x_max)
        fst_ax.set_ylabel(r'F$_{st}$', fontsize=8)
        
        # Set custom x-axis ticks only at min and max values
        fst_ax.set_xticks([x_min, x_max])
        fst_ax.set_xticklabels([f'{int(x_min):,}', f'{int(x_max):,}'])
        fst_ax.tick_params(axis='both', labelsize=6)
        fst_ax.spines[["top", "right"]].set_visible(False)

        # Add legend
        if legend_handles:
            fst_ax.legend(legend_handles, [h.get_label() for h in legend_handles],
                         bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                         mode='expand', borderaxespad=0., ncol=3, fontsize=6, frameon=False)
        
        return legend_handles

    
    def plot(self):
        """Create the gene visualization plot."""
        # Add data from config
        self.add_data(self.config['data'])
        
        # Initialize figure
        self.fig, self.ax = plt.subplots(figsize=(self.width, self.height))
        
        # Extract gene data for all genes
        for gene_id in self.genes:
            self.gene_data[gene_id] = self._extract_gene_features(gene_id)
        
        # Plot each gene according to its configuration
        for i, gene_id in enumerate(self.genes):
            gene_config = self._get_gene_config(gene_id)
            gene_label_size = gene_config.get('gene_label_size', 10)
            gene_data = self.gene_data[gene_id]
            
            # Get bbox from config or use default
            default_bbox = self.toml_config['default']['bbox']
            bbox = gene_config.get('bbox', default_bbox)
            
            # # If multiple genes, adjust bbox positions vertically
            # if len(self.genes) > 1:
            #     bbox_height = bbox[3] / len(self.genes)
            #     bbox = [bbox[0], bbox[1] + i * bbox_height, bbox[2], bbox_height]
            
            # Draw gene structure
            gene_coords = self._draw_gene_structure(gene_data, gene_config, bbox)
            
            # Draw annotations if available
            self._draw_annotations(gene_data, gene_config, bbox, gene_coords)
            
            # Draw haplotypes if available
            self._draw_haplotypes(gene_data, gene_config, bbox)
            
            # Draw pi plot if available
            self._draw_pi_plot(gene_data, gene_config, bbox)
            
            # Draw fst plot if available
            self._draw_fst_plot(gene_data, gene_config, bbox)
            
            # Add gene name label
            self.ax.annotate(
                gene_id,
                xy=(bbox[0] + bbox[2]/2, bbox[1]),
                xytext=(0, -10),
                textcoords='offset points',
                ha='center', va='center', fontsize=gene_label_size, weight='bold'
            )
        
        # Set axis properties
        # self.ax.set_xlim(-0.1, 1.1)
        # self.ax.set_ylim(-0.1, 1.1)
        self.ax.autoscale()
        self.ax.set_aspect('equal')
        
        # Remove axes and spines
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for spine in self.ax.spines.values():
            spine.set_visible(False)
        
        # Save figure
        if self.output:
            plt.savefig(self.output, dpi=300, bbox_inches='tight')
            logger.info(f"Gene plot saved to: {self.output} with size {self.width}x{self.height} inches")
        
        return self.fig, self.ax
