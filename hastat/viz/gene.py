# -*- coding: utf-8 -*-
# @Time    : 2025/07/25 12:58
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : gene.py

import csv
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from hastat.viz._export import save_figure
from hastat.viz import tracks
from hastat.utils.population import ld_r_squared
import matplotlib.patches as patches
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.legend import Legend
import tomli
from scipy.cluster import hierarchy
from hastat.log.logger import logger
from hastat.dataio.gff import read as read_gff
from hastat.dataio.gff import GeneAnnotation
from hastat.dataio.table import read_window_table
from hastat.dataio.table import read_assignments
from hastat.utils.haplotype import cluster_haplotypes


GENE_TOML_KEYS = {
    'bbox', 'layout', 'gene_y', 'gene_height', 'gene_color', 'gene_feature',
    'gene_label_size', 'gene_label_show',
    'ann_file', 'ann_y', 'ann_band_y', 'ann_legend_y', 'ann_text_size',
    'hap_file', 'hap_y', 'hap_height', 'hap_cols', 'hap_label_show', 'hap_label_size',
    'hap_count_show', 'hap_color_file', 'hap_group_bar_x', 'hap_group_bar_width',
    'hap_cluster_show', 'hap_cluster_linkage', 'hap_cluster_width',
    'hap_het_rate', 'hap_mis_rate', 'var_het_rate', 'var_mis_rate',
    'af_show', 'af_height', 'af_y', 'af_style', 'pop_file',
    'vcf_file', 'ld_show', 'ld_y', 'ld_cmap', 'ld_colorbar', 'shade_gene',
    'pi_file', 'pi_y', 'pi_height', 'pi_style',
    'fst_file', 'fst_y', 'fst_height', 'fst_style',
    'tajima_d_file', 'tajima_d_y', 'tajima_d_height', 'tajima_d_style',
}


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
        self.gene_id = None
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
        self.gene_data = None
        
        # Initialize figure and axes
        self.fig = None
        self.ax = None
        self._count_label_layout = []
        
    def add_data(self, data: dict):
        """
        Add data from the config data section.
        
        :param data: Data dictionary containing gff_file, gene, toml_file
        """
        self.gff_file = data['gff_file']
        self.gene_id = data['gene']
        if not isinstance(self.gene_id, str) or not self.gene_id.strip():
            raise ValueError('gene must be one nonempty GFF gene ID')
        self.toml_file = data['toml_file']
        
        # Load TOML configuration
        self.toml_config = self._load_toml_config()
        
        # Load GFF database
        self.gff_db = read_gff(self.gff_file)
        self.gene_feature = GeneAnnotation(self.gff_db)
        
    def _load_toml_config(self):
        """Load TOML configuration file."""
        try:
            with open(self.toml_file, 'rb') as f:
                config = tomli.load(f)
            for section, options in list(config.items()):
                if not isinstance(options, dict):
                    logger.warning("Unknown gene TOML parameter '%s' at the top level", section)
                    del config[section]
                    continue
                for key in sorted(options.keys() - GENE_TOML_KEYS):
                    logger.warning("Unknown gene TOML parameter '%s' in [%s]", key, section)
                    del options[key]
            logger.info(f"Successfully loaded TOML configuration: {self.toml_file}")
            return config
        except Exception as e:
            logger.error(f"Failed to load TOML configuration: {e}")
            raise e
    
    def _get_gene_config(self, gene_id):
        """Get configuration for a specific gene, falling back to default if not found."""
        defaults = self.toml_config.get('default', {})
        selected = self.toml_config.get(gene_id, {})
        if not isinstance(defaults, dict) or not isinstance(selected, dict):
            raise ValueError('Gene TOML default and gene-specific entries must be tables')
        gene_config = dict(defaults, **selected)
        if 'bbox' in gene_config:
            bbox = gene_config['bbox']
            if (not isinstance(bbox, list) or len(bbox) != 4
                    or any(type(v) not in (int, float) or not np.isfinite(v) for v in bbox)
                    or min(bbox[2:]) <= 0):
                raise ValueError('bbox requires [x, y, width, height] with finite values and positive dimensions')
        mode = gene_config.get('layout', 'fixed')
        if mode not in ('fixed', 'auto'):
            raise ValueError('layout must be fixed or auto')
        if mode == 'auto':
            manual = sorted({'pi_y', 'fst_y', 'tajima_d_y', 'af_y', 'ld_y'} & gene_config.keys())
            if manual:
                raise ValueError('Use layout="fixed" for explicit positions: ' + ', '.join(manual))
        for key in ('gene_y', 'gene_height', 'ann_band_y', 'ann_legend_y', 'ld_y'):
            if key in gene_config:
                value = gene_config[key]
                if type(value) not in (int, float) or not np.isfinite(value):
                    raise ValueError(f'{key} must be a finite bbox-relative coordinate')
        for option in ('ld_show', 'ld_colorbar', 'shade_gene', 'af_show'):
            if option in gene_config and not isinstance(gene_config[option], bool):
                raise ValueError(f'{option} must be true or false')
        if gene_config.get('ld_show', False) and not gene_config.get('vcf_file'):
            raise ValueError('ld_show=true requires an indexed vcf_file')
        if gene_config.get('hap_file'):
            prefix = Path(gene_config['hap_file']).expanduser()
            if not prefix.is_absolute():
                prefix = Path(self.toml_file).resolve().parent / prefix
            # Also accept an explicit .table.csv path; there is still only one input setting.
            prefix = str(prefix)
            if prefix.endswith('.table.csv'):
                prefix = prefix[:-len('.table.csv')]
            gene_config['hap_file'] = prefix + '.table.csv'
            group_path = Path(prefix + '.group.csv')
            if group_path.exists() or group_path.is_symlink():
                if not group_path.is_file():
                    raise ValueError(f'Haplotype sample groups must be a readable file: {group_path}')
                gene_config['_hap_sample_file'] = str(group_path)
        for name in ('ann_file', 'hap_file', 'hap_color_file', 'pop_file', 'vcf_file', 'pi_file', 'fst_file', 'tajima_d_file'):
            value = gene_config.get(name)
            if value:
                path = Path(value).expanduser()
                if not path.is_absolute():
                    path = Path(self.toml_file).resolve().parent / path
                if not path.is_file():
                    raise FileNotFoundError(f'Configured {name} does not exist: {path}')
                gene_config[name] = str(path)
        return gene_config
    
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
        return tracks.read_annotations(ann_file)

    def _load_haplotype_data(self, hap_file, hap_cols=None):
        """Load haplotype data from CSV file."""
        if not hap_file or hap_file == "":
            return None
        try:
            with open(hap_file, encoding='utf-8-sig', newline='') as stream:
                header = next(csv.reader(stream), [])
            meta = ['chrom', 'pos', 'ref', 'alt']
            if not set(meta).issubset(header) or len(header) != len(set(header)) or any(not name.strip() for name in header):
                raise ValueError('Haplotype table requires chrom,pos,ref,alt and unique nonempty column names')
            hap_data = pd.read_csv(hap_file, dtype=str, keep_default_na=False)
            if hap_data[meta].apply(lambda col: col.str.strip().eq('')).any().any():
                raise ValueError('Haplotype variant identifiers must be nonempty')
            positions = pd.to_numeric(hap_data['pos'], errors='raise')
            if not (np.isfinite(positions) & (positions >= 1) & (positions == np.floor(positions))).all():
                raise ValueError('Haplotype positions must be positive integers')
            hap_data['pos'] = positions.astype('int64')
            if hap_data.duplicated(meta).any():
                raise ValueError('Duplicate variant records in haplotype table')
            all_haps = [col for col in header if col not in meta]
            for hap in all_haps:
                values = hap_data[hap].replace({'': '-1', 'NA': '-1', 'NaN': '-1', 'nan': '-1'})
                values = pd.to_numeric(values, errors='raise')
                if not values.isin([-1, 0, 1, 2]).all():
                    raise ValueError('Haplotype genotypes must be -1 (missing), 0, 1 or 2')
                hap_data[hap] = values.astype('int8')
            logger.info(f"Loaded haplotype data from {hap_file}")
            if isinstance(hap_cols, list) and len(hap_cols) > 0:
                hap_data = hap_data.loc[:, ['chrom', 'pos', 'ref', 'alt'] + hap_cols]
            hap_data.attrs['all_haplotypes'] = all_haps
            logger.info(f"detected {len(hap_data.columns) - 4} haplotypes with {len(hap_data)} variant sites")
            return hap_data
        except Exception as e:
            raise ValueError(f"Failed to load haplotype file {hap_file}: {e}") from e
    
    def _load_haplotype_group_data(self, group_file):
        """Load haplotype group data from CSV file."""
        if not group_file or group_file == "":
            return None
        try:
            group_data = pd.read_csv(group_file, header=None, dtype=str, keep_default_na=False)
            # Assuming 2 columns: haplotype, group
            if len(group_data.columns) == 2:
                group_data.columns = ['haplotype', 'group']
                group_data = group_data.apply(lambda col: col.str.strip())
                if group_data.eq('').any().any() or group_data.haplotype.duplicated().any():
                    raise ValueError('hap_color_file requires one nonempty category per unique haplotype; duplicate ID found or empty value')
                logger.info(f"Loaded haplotype group data from {group_file}")
                return group_data
            else:
                raise ValueError(f"Invalid group file format: {group_file}")
        except Exception as e:
            raise ValueError(f"Failed to load group file {group_file}: {e}") from e
    
    def _physical_frame(self, gene_config, bbox):
        """Horizontal frame shared by all tracks using physical genomic positions."""
        x, y, width, height = bbox
        return x, y + gene_config.get('gene_y', 0.5) * height, width, gene_config.get('gene_height', 0.2) * height

    def _coordinate_transform(self, pos, bbox, region_start, region_end):
        """Map base centers/edges on a 1-based inclusive region, including one base."""
        x, _, width, _ = bbox
        return x + (pos - (region_start - 0.5)) / (region_end - region_start + 1) * width
    
    def _draw_gene_structure(self, gene_data, gene_config, bbox):
        """Draw gene structure features (CDS, exons, UTRs)."""
        x_start, y_start, bbox_width, bbox_height = bbox
        color = gene_config.get('gene_color', '#24495D')
        feature_type = gene_config.get('gene_feature', 'CDS')
        gene_label_size = gene_config.get('gene_label_size', 10)
        gene_label_show = gene_config.get('gene_label_show', True)
        
        # Calculate absolute coordinates
        frame = self._physical_frame(gene_config, bbox)
        abs_x, abs_y, abs_width, abs_height = frame
        
        gene_start = gene_data['gene_start']
        gene_end = gene_data['gene_end']
        region_start, region_end = gene_data['start'], gene_data['end']
        gene_left = self._coordinate_transform(gene_start - 0.5, frame, region_start, region_end)
        gene_right = self._coordinate_transform(gene_end + 0.5, frame, region_start, region_end)
        
        # Draw selected features
        features_to_draw = gene_data.get(feature_type, [])
        if len(features_to_draw) == 0:
            logger.warning(f"No {feature_type} features found! Please check the GFF file and toml file.")
        for feature_start, feature_end in features_to_draw:
            # Inclusive features occupy their full base cells, even if start == end.
            feature_x = self._coordinate_transform(feature_start - 0.5, frame, region_start, region_end)
            feature_right = self._coordinate_transform(feature_end + 0.5, frame, region_start, region_end)
            feature_width = feature_right - feature_x
            
            feature_rect = patches.Rectangle(
                (feature_x, abs_y - abs_height/2),
                feature_width, abs_height,
                facecolor=color, edgecolor=color, 
                alpha=1, zorder=2
            )
            self.ax.add_patch(feature_rect)
        
        # Connect features behind their fills and mark direction once at the terminal edge.
        if features_to_draw:
            backbone = patches.PathPatch(
                mpl.path.Path([(gene_left, abs_y), (gene_right, abs_y)]),
                fill=False, edgecolor='#263746', linewidth=.9, zorder=1)
            self.ax.add_patch(backbone)
            if gene_data['strand'] in ('+', '-'):
                sign = 1 if gene_data['strand'] == '+' else -1
                terminal = (gene_right if sign == 1 else gene_left, abs_y)
                marker = self.ax.annotate('>' if sign == 1 else '<', xy=terminal,
                    xytext=(3*sign, 0), textcoords='offset points',
                    ha='left' if sign == 1 else 'right', va='center',
                    fontsize=12, color='#263746', annotation_clip=False, zorder=2)
                marker.set_gid('strand-direction')

        # Add gene name label
        if gene_label_show:
            self.ax.annotate(
                gene_data['gene_id'],
                xy=((gene_left + gene_right)/2, abs_y + abs_height/2),
                xytext=(0, 5),
                textcoords='offset points',
                ha='center', va='bottom', fontsize=gene_label_size, weight='bold',
                zorder=10
            )

        return abs_x, abs_y, abs_width, abs_height
    
    def _prepare_haplotypes(self, gene_data, gene_config):
        """Draw haplotype information."""
        hap_file = gene_config.get('hap_file', '')
        hap_cols = gene_config.get('hap_cols', [])
        sample_hap_file = gene_config.get('_hap_sample_file', '')
        count_show = gene_config.get('hap_count_show', bool(sample_hap_file))
        if not isinstance(count_show, bool):
            raise ValueError('hap_count_show must be true or false')
        assignments = read_assignments(sample_hap_file, label='haplotypes') if sample_hap_file else None
        hap_counts = assignments.haplotypes.value_counts().to_dict() if count_show and assignments is not None else None
        
        # New config options
        hap_color_file = gene_config.get('hap_color_file', '')
        # Clustering config
        hap_cluster_show = gene_config.get('hap_cluster_show', False)
        hap_cluster_linkage = gene_config.get('hap_cluster_linkage', 'average')
        
        hap_data = self._load_haplotype_data(hap_file)
        all_profiles = hap_data
        if hap_data is None:
            return
        if assignments is not None:
            unknown = sorted(set(assignments.haplotypes) - set(hap_data.attrs['all_haplotypes']))
            if unknown:
                raise ValueError('Sample assignments contain haplotypes absent from the table: ' + ', '.join(unknown))

        # Filter haplotypes and variants
        hap_het_rate = gene_config.get('hap_het_rate', [0, 1])
        hap_mis_rate = gene_config.get('hap_mis_rate', [0, 1])
        var_het_rate = gene_config.get('var_het_rate', [0, 1])
        var_mis_rate = gene_config.get('var_mis_rate', [0, 1])

        # 1. Filter haplotypes (columns)
        meta_cols = ['chrom', 'pos', 'ref', 'alt']
        current_hap_cols = hap_cols or [c for c in hap_data.columns if c not in meta_cols]
        if not set(current_hap_cols) <= set(hap_data.columns) - set(meta_cols):
            raise ValueError('Requested haplotype absent from the table')
        
        valid_haps = []
        for hap in current_hap_cols:
            series = hap_data[hap]
            n_total = len(series)
            n_miss = (series == -1).sum()
            n_called = n_total - n_miss
            
            mis_rate = n_miss / n_total if n_total > 0 else 0
            het_rate = (series == 1).sum() / n_called if n_called > 0 else 0
            
            if (hap_mis_rate[0] <= mis_rate <= hap_mis_rate[1]) and \
               (hap_het_rate[0] <= het_rate <= hap_het_rate[1]):
                valid_haps.append(hap)
        
        if not valid_haps:
            logger.warning("No haplotypes left after filtering!")
            return

        hap_data = hap_data[meta_cols + valid_haps]
        current_hap_cols = valid_haps

        # 2. Filter variants (rows)
        haps_df = hap_data[current_hap_cols]
        n_haps = len(current_hap_cols)
        
        if n_haps > 0:
            n_miss_per_site = (haps_df == -1).sum(axis=1)
            n_called_per_site = n_haps - n_miss_per_site
            
            mis_rate_per_site = n_miss_per_site / n_haps
            
            het_counts = (haps_df == 1).sum(axis=1)
            het_rate_per_site = pd.Series(0.0, index=hap_data.index)
            mask_called = n_called_per_site > 0
            het_rate_per_site[mask_called] = het_counts[mask_called] / n_called_per_site[mask_called]
            
            mask_keep = (mis_rate_per_site >= var_mis_rate[0]) & (mis_rate_per_site <= var_mis_rate[1]) & \
                        (het_rate_per_site >= var_het_rate[0]) & (het_rate_per_site <= var_het_rate[1])
            
            hap_data = hap_data[mask_keep]
            
            if hap_data.empty:
                logger.warning("No variants left after filtering!")
                return
            
            logger.info(f"After filtering: {len(hap_data.columns) - 4} haplotypes, {len(hap_data)} sites")
        
        
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

        # Load group data
        group_data = self._load_haplotype_group_data(hap_color_file)
        group_info = {} # hap -> group
        if group_data is not None:
            for _, row in group_data.iterrows():
                group_info[row['haplotype']] = str(row['group'])
        
        # Sort haplotypes
        if not hap_cluster_show:
            if group_info:
                # Filter haplotypes that are in group_info
                hap_cols = [hap for hap in hap_cols if hap in group_info]
                
                if not hap_cols:
                    logger.warning("No haplotypes found in group file!")
                    return

                # Sort by group then by haplotype name
                hap_cols.sort(key=lambda x: (group_info[x], x))
            else:
                hap_cols.sort()
        else:
            if group_info:
                hap_cols = [hap for hap in hap_cols if hap in group_info]

        if not hap_cols:
            logger.warning('No haplotypes left to display after group selection')
            return
        if assignments is not None:
            assigned_haps = set(assignments.haplotypes)
            missing = [hap for hap in hap_cols if hap not in assigned_haps]
            if missing:
                raise ValueError('No sample assignments for displayed haplotypes: ' + ', '.join(missing))
        hap_labels = {hap: f'{hap} (n={int(hap_counts[hap])})' if hap_counts is not None else hap
                      for hap in hap_cols}

        # Sort gene_hap by position to ensure consistent column order
        gene_hap_sorted = gene_hap.sort_values(['pos','ref','alt'])

        # Perform clustering if enabled
        Z = None
        if hap_cluster_show:
            hap_cols, Z = cluster_haplotypes(gene_hap_sorted, hap_cols, hap_cluster_linkage)

        keys = pd.MultiIndex.from_frame(gene_hap_sorted[tracks.KEY])
        sites = all_profiles.set_index(tracks.KEY).loc[keys].reset_index()
        genotypes = None
        show_ld = gene_config.get('ld_show', False)
        show_af = gene_config.get('af_show', True)
        if gene_config.get('vcf_file') and (show_ld or show_af):
            genotypes = tracks.aligned_genotypes(gene_config['vcf_file'], sites)
        freq = None
        if show_af:
            freq = tracks.allele_frequencies(sites, assignments, gene_config.get('pop_file'), genotypes)
        ld = ld_r_squared(genotypes.to_numpy()) if show_ld and len(sites)>1 else None
        return dict(sites=sites, haps=hap_cols, labels=hap_labels, groups=group_info, linkage=Z,
                    frequencies=freq, ld=ld)

    def _draw_haplotypes(self, gene_data, gene_config, bbox, view, freq_bottom, ld_top):
        """Draw existing haplotype geometry and the aligned data tracks."""
        gene_hap_sorted = view['sites']
        hap_cols, hap_labels, group_info, Z = view['haps'], view['labels'], view['groups'], view['linkage']
        x_start, y_start, bbox_width, bbox_height = bbox
        hap_y = gene_config.get('hap_y', 0)
        hap_height = gene_config.get('hap_height', 0.3)
        hap_label_show = gene_config.get('hap_label_show', True)
        hap_label_size = gene_config.get('hap_label_size', 8)
        af_height = gene_config.get('af_height', 0.05)
        af_style = gene_config.get('af_style', 'bar')
        hap_group_bar_x = gene_config.get('hap_group_bar_x', -0.12)
        hap_group_bar_width = gene_config.get('hap_group_bar_width', 0.02)
        
        hap_cluster_width = gene_config.get('hap_cluster_width', 0.1)
        # Calculate coordinates for equal-spaced grid
        hap_abs_x = x_start
        hap_abs_y = y_start + hap_y * bbox_height
        hap_abs_width = bbox_width
        hap_abs_height = hap_height * bbox_height
        
        num_haps = len(hap_cols)  # Number of rows (haplotypes)
        num_sites = len(gene_hap_sorted)
        
        # Calculate cell dimensions
        cell_width = hap_abs_width / num_sites
        cell_height = hap_abs_height / num_haps
        
        # Draw haplotype matrix with equal spacing
        colors = {0: 'white', 1: '#90CAEE', 2: '#114F8B', -1: 'lightgray'}
        
        # Determine edge properties based on number of haplotypes
        edge_color = 'black' if num_haps <= 10 else None
        line_width = 0.5 if num_haps <= 10 else 0
        
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
                    facecolor=color, edgecolor=edge_color, linewidth=line_width
                )
                self.ax.add_patch(rect)
        
        # Draw labels and group bars
        label_artists, group_bar_artists, dendrogram_artists = [], [], []
        if group_info:
            # Group colors
            unique_groups = sorted(list(set(group_info.values())))
            # Use a colormap
            group_cmap = LinearSegmentedColormap.from_list('group_cmap', ['#C5504B', '#114F8B', '#FCE988', '#90CAEE'], N=100)
            colors = group_cmap(np.linspace(0, 1, len(unique_groups)))
            group_colors = {grp: colors[i] for i, grp in enumerate(unique_groups)}

            for i, hap in enumerate(hap_cols):
                grp = group_info[hap]
                color = group_colors[grp]
                
                # Draw group color block for this row
                rect_x = hap_abs_x + hap_group_bar_x
                rect_y = hap_abs_y + i * cell_height
                
                rect = patches.Rectangle(
                    (rect_x, rect_y),
                    hap_group_bar_width, cell_height,
                    facecolor=color, edgecolor=None, linewidth=0
                )
                self.ax.add_patch(rect)
                group_bar_artists.append(rect)
                
                # Draw hap label
                if hap_label_show:
                    label_y = hap_abs_y + i * cell_height + cell_height/2
                    label_artists.append(self.ax.text(
                        hap_abs_x - 0.02, label_y,
                        hap_labels[hap], ha='right', va='center', fontsize=hap_label_size
                    ))

            displayed = sorted({group_info[hap] for hap in hap_cols})
            handles = [patches.Patch(facecolor=group_colors[group]) for group in displayed]
            # Keep the legend below the matrix and to the left of the AF and LD tracks.
            top = hap_abs_y
            if view['frequencies'] is not None:
                top = min(top, freq_bottom)
            if view['ld'] is not None:
                top = min(top, ld_top)
            right = hap_abs_x + min(0, hap_group_bar_x + hap_group_bar_width)
            legend = Legend(self.ax, handles, displayed, loc='upper right',
                            bbox_to_anchor=(right, top), bbox_transform=self.ax.transData,
                            frameon=False, fontsize=8, handlelength=1, handletextpad=.5,
                            borderpad=0, borderaxespad=.8, labelspacing=.4)
            legend.set_gid('hap-color-legend')
            self.ax.add_artist(legend)
        else:
            # Add haplotype labels on the left
            if hap_label_show:
                for hap_idx, hap_col in enumerate(hap_cols):
                    label_y = hap_abs_y + hap_idx * cell_height + cell_height/2
                    label_artists.append(self.ax.text(
                        hap_abs_x - 0.02, label_y,
                        hap_labels[hap_col], ha='right', va='center', fontsize=hap_label_size
                    ))
        
        # Draw dendrograms
        if Z is not None:
            # Determine dendrogram position
            dendro_x_end = hap_abs_x + hap_group_bar_x - 0.02 if group_info else hap_abs_x - 0.05
            dendro_width = hap_cluster_width
            
            dendro = hierarchy.dendrogram(Z, no_plot=True)
            icoord = np.array(dendro['icoord'])
            dcoord = np.array(dendro['dcoord'])
            
            # Map icoord (leaves) to y
            # scipy dendrogram leaves are at 5, 15, 25...
            # We map them to the center of our rows
            # normalized_y = (y - 5) / 10  -> 0, 1, 2...
            
            y_coords = hap_abs_y + ((icoord - 5)/10) * cell_height + cell_height/2
            
            # Map dcoord to x
            max_d = np.max(dcoord) if np.max(dcoord) > 0 else 1
            x_coords = dendro_x_end - (dcoord / max_d) * dendro_width
            
            # Plot lines
            for xs, ys in zip(x_coords, y_coords):
                dendrogram_artists.extend(self.ax.plot(xs, ys, color='black', linewidth=0.5))

        if label_artists:
            self._count_label_layout.append((label_artists, group_bar_artists, dendrogram_artists))

        freq = view['frequencies']
        if freq is not None:
            tracks.draw_frequencies(self.ax, freq, hap_abs_x, freq_bottom,
                                         bbox_width, af_height*bbox_height, af_style)
        if view['ld'] is not None:
            tracks.draw_ld(self.ax, view['ld'], hap_abs_x, ld_top, bbox_width,
                                gene_config.get('ld_cmap', 'Blues'), gene_config.get('ld_colorbar', True))

        # Draw connecting lines from haplotype columns to gene variant positions
        self._draw_haplotype_connections(gene_data, gene_hap_sorted, self._physical_frame(gene_config, bbox),
                                       hap_abs_x, hap_abs_y, hap_abs_height, 
                                       cell_width)
    
    def _fit_haplotype_count_labels(self):
        """Keep group bars/trees left of the measured count labels without shrinking the matrix."""
        if not self._count_label_layout:
            return
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        padding = self.fig.dpi * 4 / 72  # four points, independent of export resolution
        for labels, bars, tree in self._count_label_layout:
            decorations = bars if bars else tree
            if not decorations:
                continue
            label_left = min(label.get_window_extent(renderer).x0 for label in labels)
            right = max(artist.get_window_extent(renderer).x1 for artist in decorations)
            shift_pixels = min(0., label_left - padding - right)
            if shift_pixels == 0:
                continue
            inverse = self.ax.transData.inverted()
            shift = inverse.transform((shift_pixels, 0))[0] - inverse.transform((0, 0))[0]
            for bar in bars:
                bar.set_x(bar.get_x() + shift)
                bar.set_clip_on(False)
            for line in tree:
                line.set_xdata(np.asarray(line.get_xdata()) + shift)
                line.set_clip_on(False)

    def _draw_haplotype_connections(self, gene_data, gene_hap_sorted, frame,
                                  hap_abs_x, hap_abs_y, hap_abs_height, cell_width):
        """Draw connecting lines from haplotype columns to gene variant positions."""
        _, gene_abs_y, _, gene_abs_height = frame
        
        # Draw connection lines for each variant site
        for site_idx, (_, row) in enumerate(gene_hap_sorted.iterrows()):
            pos = row['pos']
            
            # Position in haplotype matrix (top of the column)
            hap_column_x = hap_abs_x + site_idx * cell_width + cell_width/2
            hap_column_y_top = hap_abs_y + hap_abs_height  # Top of haplotype matrix
            
            # Position on gene (variant location)
            gene_variant_x = self._coordinate_transform(pos, frame, gene_data['start'], gene_data['end'])
            gene_variant_y = gene_abs_y - gene_abs_height/2  # Bottom edge of features
            
            # Draw connecting line: diagonal line first, then vertical line to gene bottom
            # Connection height is half the distance between haplotype top and gene bottom
            intermediate_y = (gene_variant_y + hap_column_y_top) / 2

            # Diagonal line from haplotype column to intermediate point
            self.ax.plot([hap_column_x, gene_variant_x], 
                        [hap_column_y_top, intermediate_y],
                        color='gray', alpha=0.6, linewidth=1, linestyle='-')
            
            # Vertical line up to gene variant position (bottom of gene)
            self.ax.plot([gene_variant_x, gene_variant_x], 
                        [intermediate_y, gene_variant_y],
                        color='gray', alpha=0.6, linewidth=1, linestyle='-')

    def _load_pi_data(self, pi_file):
        return read_window_table(pi_file, 'pi') if pi_file else None

    def _load_fst_data(self, fst_file):
        return read_window_table(fst_file, 'fst') if fst_file else None

    def _load_tajima_d_data(self, path):
        return read_window_table(path, 'tajima_d') if path else None

    def _draw_pi_plot(self, gene_data, gene_config, bbox):
        return self._draw_window_track(gene_data, gene_config, bbox, 'pi')

    def _draw_fst_plot(self, gene_data, gene_config, bbox):
        return self._draw_window_track(gene_data, gene_config, bbox, 'fst')

    def _draw_tajima_d_plot(self, gene_data, gene_config, bbox):
        return self._draw_window_track(gene_data, gene_config, bbox, 'tajima_d')

    def _draw_window_track(self, gene_data, gene_config, bbox, statistic):
        """Draw existing window estimates with literal labels and visible points."""
        path = gene_config.get(f'{statistic}_file', '')
        frame = read_window_table(path, statistic) if path else None
        if frame is None or frame.empty:
            return None
        style = gene_config.get(f'{statistic}_style', 'fill')
        if style not in ('line', 'fill'):
            raise ValueError(f'{statistic}_style must be line or fill')
        if 'chrom' in frame:
            frame = frame.loc[frame.chrom.eq(str(gene_data['chrom']))].copy()
            if frame.empty:
                logger.warning(f'No {statistic} windows on chromosome {gene_data["chrom"]}')
                return None
        start, end = int(frame.start.min()), int(frame.end.max())
        frame['pos_bp'] = frame.start + (frame.end - frame.start) / 2
        x, y, width, height = bbox
        track_y = gene_config.get(f'{statistic}_y', {'pi': 0.55, 'fst': 0.8, 'tajima_d': 0.3}[statistic])
        track_height = gene_config.get(f'{statistic}_height', 0.2)
        track_ax = self.ax.inset_axes([x, y + track_y * height, width, track_height * height],
                                     transform=self.ax.transData)
        if gene_config.get('shade_gene', False):
            left, right = max(start,gene_data['gene_start']), min(end,gene_data['gene_end'])
            if left <= right:
                # Genomic positions are 1-based inclusive; shade the complete base cells.
                track_ax.axvspan(left-.5,right+.5,facecolor='#E8E8E8',edgecolor='none',
                                 zorder=-2, gid='gene-span')
        cmap = LinearSegmentedColormap.from_list(f'{statistic}_cmap',
            ['#C5504B', '#114F8B', '#FCE988', '#90CAEE'], N=100)
        axis_label = {'pi': r'$\mathit{\pi}$ ($10^{-3}$)', 'fst': r'$F_{\mathrm{ST}}$', 'tajima_d': "Tajima's $D$"}[statistic]
        if statistic == 'tajima_d':
            track_ax.axhline(0, color='0.6', linewidth=0.6, linestyle='--', zorder=0)
        if 'groups' in frame:
            labels = sorted(frame.groups.unique())
            colors = cmap(np.linspace(0, 1, len(labels)))
            series = [(label, frame.loc[frame.groups.eq(label)], color)
                      for label, color in zip(labels, colors)]
        else:
            series = [(axis_label, frame,
                       'C0' if statistic == 'pi' else 'C1')]
        handles = []
        for label, group, color in series:
            group = group.sort_values('pos_bp', kind='stable')
            positions = group.pos_bp.to_numpy(dtype=float)
            values = group[statistic].to_numpy(dtype=float) * (1000 if statistic == 'pi' else 1)
            if not np.isfinite(values).any():
                logger.warning(f'No finite {statistic} estimates for group {label} in this region')
            if style == 'fill':
                track_ax.fill_between(positions, 0, values, color=color, alpha=0.6)
            handles.append(track_ax.plot(positions, values, color=color, linewidth=1,
                marker='o', markersize=2.5, label=label)[0])
        # Even a single-base window occupies its full inclusive base cell.
        track_ax.set_xlim(start - .5, end + .5)
        ticks = [start, end] if start < end else [start]
        track_ax.set_xticks(ticks)
        track_ax.set_xticklabels([f'{value:,}' for value in ticks])
        track_ax.set_ylabel(axis_label, fontsize=8)
        track_ax.tick_params(axis='both', labelsize=6)
        track_ax.spines[['top', 'right']].set_visible(False)
        track_ax.legend(handles, [handle.get_label() for handle in handles],
            bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', mode='expand',
            borderaxespad=0., ncol=3, fontsize=6, frameon=False)
        return handles

    def _draw_group_bar(self, group, start_idx, end_idx, hap_abs_x, hap_abs_y, cell_height, color, x_offset, width):
        """Draw group bar and label."""
        # y range: from start_idx * cell_height to (end_idx + 1) * cell_height
        y_bottom = hap_abs_y + start_idx * cell_height
        height = (end_idx - start_idx + 1) * cell_height
        
        # Position: left of hap label.
        # Assuming hap label is around hap_abs_x - 0.02 to -0.08 depending on length.
        # Let's put Group Rect at hap_abs_x + x_offset
        rect_width = width
        rect_x = hap_abs_x + x_offset
        
        rect = patches.Rectangle(
            (rect_x, y_bottom),
            rect_width, height,
            facecolor=color, edgecolor='black', linewidth=0.5
        )
        self.ax.add_patch(rect)
        
        # Group Label
        label_y = y_bottom + height / 2
        self.ax.text(
            rect_x - 0.01, label_y,
            group, ha='right', va='center', fontsize=8, fontweight='bold'
        )
    
    @staticmethod
    def _track_layout(conf, bbox, view):
        """Keep configured positions in bbox coordinates and attach LD below AF or the matrix."""
        _, y, width, height = bbox
        freq_y = y + conf.get('af_y', 0) * height
        layout = dict(gene_y=y, freq_y=freq_y, ld_top=y, height=height, windows={})
        has_ld = view is not None and view['ld'] is not None
        freq = view['frequencies'] if view is not None else None
        freq_height = 0 if freq is None else len(freq.columns)*conf.get('af_height',.05)*height
        if conf.get('layout', 'fixed') == 'fixed':
            if not has_ld:
                return layout
            matrix_bottom = y + conf.get('hap_y', 0) * height
            if freq is not None and 'af_y' not in conf:
                freq_y = matrix_bottom - .012 * width - freq_height
            layout['freq_y'] = freq_y
            layout['ld_top'] = (freq_y if freq is not None else matrix_bottom) - .014 * width
            if 'ld_y' in conf:
                layout['ld_top'] = y + conf['ld_y'] * height
            return layout
        ld_top = y + width/2 if has_ld else y
        freq_y = ld_top + .014*width if has_ld else y
        matrix_bottom = freq_y+freq_height+.012*width
        gene_y = matrix_bottom-conf.get('hap_y',0)*height if view is not None else y
        gene_top = gene_y+(conf.get('gene_y',.5)+conf.get('gene_height',.2)/2)*height
        legend_gap = .045*width if conf.get('ann_file') else 0.
        labels_top = max(gene_top+.1*width, gene_y+conf.get('ann_y',0)*height+.04*width) + legend_gap
        cursor = labels_top
        windows = {}
        for stat in ('tajima_d','fst','pi'):
            if conf.get(stat+'_file'):
                bottom = cursor+.04*width  # space for physical-coordinate tick labels
                windows[stat+'_y'] = (bottom-y)/height
                cursor = bottom+conf.get(stat+'_height',.2)*height+.04*width  # series legend
        layout.update(gene_y=gene_y,freq_y=freq_y,ld_top=ld_top,
                      height=max(height,cursor-y+.01*width),windows=windows)
        return layout

    def plot(self):
        """Create the gene visualization plot."""
        self._count_label_layout = []
        # Add data from config
        self.add_data(self.config['data'])
        
        # Initialize figure
        self.fig, self.ax = plt.subplots(figsize=(self.width, self.height))
        
        conf = self._get_gene_config(self.gene_id)
        self.gene_data = self._extract_gene_features(self.gene_id)
        bbox = conf.get('bbox', [0.1, 0.1, 0.8, 0.8])
        view = self._prepare_haplotypes(self.gene_data, conf)
        if conf.get('ld_show', False) and view is None:
            raise ValueError(f'LD requires a haplotype table with displayed sites for {self.gene_id}')
        layout = self._track_layout(conf, bbox, view)
        delta = self._render_gene(conf, bbox, view, layout)
        if delta:
            for key in layout['windows']:
                layout['windows'][key] -= delta / bbox[3]
            layout['height'] = max(bbox[3], layout['height'] - delta)
            self.ax.clear()
            self._count_label_layout = []
            self.fig.set_figheight(self.height)
            self._render_gene(conf, bbox, view, layout, measure=False)

        if self.output:
            save_figure(self.fig, self.output, dpi=300, bbox_inches='tight')
            logger.info(f'Gene plot saved to: {self.output}')
        return self.fig, self.ax

    def _render_gene(self, conf, bbox, view, layout, measure=True):
        gene = self.gene_data
        x, y, w, h = bbox
        upper = (x, layout['gene_y'], w, h)
        frame = self._draw_gene_structure(gene, conf, upper)
        grid = None
        if view is not None:
            self._draw_haplotypes(gene, conf, upper, view, layout['freq_y'], layout['ld_top'])
            band_y = upper[1]+(conf.get('hap_y',0)+conf.get('hap_height',.3))*h+.004*w
            grid = (x, band_y, w)
        transform = lambda pos: self._coordinate_transform(pos, frame, gene['start'], gene['end'])
        ann = self._load_annotation_data(conf.get('ann_file'))
        anchors, legend = [], None
        if ann is not None:
            anchors, legend = tracks.draw_annotations(self.ax, ann, gene, frame,
                view['sites'] if view else None, grid, transform,
                band_bottom=y+conf['ann_band_y']*h if 'ann_band_y' in conf else None)
        label_y = upper[1]+conf['ann_y']*h if 'ann_y' in conf else None
        title = next((t for t in reversed(self.ax.texts) if t.get_text()==gene['gene_id']), None)
        center = transform((gene['gene_start']+gene['gene_end'])/2)
        anchor = (center, frame[1]+frame[3]/2)
        window_config = dict(conf, **layout['windows'])
        windows = []
        for stat in ('pi', 'fst', 'tajima_d'):
            lines = self._draw_window_track(gene, window_config, bbox, stat)
            if lines and layout['windows']:
                windows.append(lines[0].axes)
        legend_y = y+conf['ann_legend_y']*h if 'ann_legend_y' in conf else None
        self.ax.update_datalim([(x,y),(x+w,y+layout['height'])])
        self.ax.autoscale()
        growth = layout['height'] - h
        if growth:
            added = growth * (1 + 2*self.ax.margins()[1])
            bounds = self.ax.get_position(original=True)
            xspan = np.ptp(self.ax.get_xlim())
            yspan = np.ptp(self.ax.get_ylim()) - added
            inches_per_unit = min(self.width*bounds.width/xspan,
                                  self.height*bounds.height/max(yspan,1e-9))
            self.fig.set_figheight(self.height+added*inches_per_unit/bounds.height)
        self.ax.set_aspect('equal')
        
        # Remove axes and spines
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for spine in self.ax.spines.values():
            spine.set_visible(False)

        self.fig.canvas.draw()
        delta = 0.
        renderer = self.fig.canvas.get_renderer()
        labels = tracks.place_labels(self.ax, anchors, conf.get('ann_text_size',8), label_y)
        if legend is not None:
            if legend_y is None:
                tracks.position_legend(self.ax, legend, anchor, [*labels,title])
            else:
                legend.set_bbox_to_anchor((anchor[0], legend_y), transform=self.ax.transData)
        if measure and windows:
            # Leave 12 points between the highest label and the lowest window axis.
            top = self.ax.transData.transform(anchor)[1]
            for artist in [*labels,title,legend]:
                if artist is not None:
                    top = max(top, artist.get_window_extent(renderer).ymax)
            bottom = min(axis.get_tightbbox(renderer).ymin for axis in windows)
            delta_px = bottom - top - 12 * self.fig.dpi / 72
            if abs(delta_px) > .5:
                inverse = self.ax.transData.inverted()
                delta = inverse.transform((0,delta_px))[1] - inverse.transform((0,0))[1]
        self._fit_haplotype_count_labels()
        return delta
