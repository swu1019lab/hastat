# -*- coding: utf-8 -*-
# @Time    : 2024/9/30 18:07
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : plot.py

import time
from hastat.log.logger import logger
from hastat.viz import bar, pie, box, gene2, network
from hastat.manage import Manager
import matplotlib as mpl
mpl.use('Agg')

# Initialize the manager and register plot types
manager = Manager()
manager.register('bar', bar.HapBar)
manager.register('pie', pie.HapPie)
manager.register('box', box.HapBox)
manager.register('gene', gene2.FancyGene)
manager.register('network', network.FancyNetwork)


def run(args):
    logger.info('Start to run plot function')
    start = time.time()

    # Get plot type from args
    plot_type = getattr(args, 'plot_type', None)
    if plot_type is None:
        logger.error('No plot type specified')
        return

    # Create config dictionary from command line arguments based on plot type
    if plot_type == 'bar':
        config = {
            'data': {
                'sample_hap': args.sample_hap,
                'sample_group': args.sample_group
            },
            'plot': {
                'group_index': args.group_index,
                'calc_percentage': args.calc_percentage,
                'haplotypes': args.haplotypes or [],
                'x_label': args.x_label,
                'y_label': args.y_label,
                'save_fig': args.output,
                'width': args.width,
                'height': args.height,
            }
        }
    elif plot_type == 'pie':
        config = {
            'data': {
                'sample_hap': args.sample_hap,
                'sample_group': args.sample_group
            },
            'plot': {
                'group_index': args.group_index,
                'calc_percentage': args.calc_percentage,
                'haplotypes': args.haplotypes or [],
                'width': args.width,
                'height': args.height,
                'save_fig': args.output,
            }
        }
    elif plot_type == 'box':
        config = {
            'data': {
                'sample_hap': args.sample_hap,
                'sample_phe': args.sample_phe
            },
            'plot': {
                'phe_index': args.phe_index,
                'width': args.width,
                'height': args.height,
                'comparisons': args.comparisons or [],
                'method': args.method,
                'step_size': args.step_size or [],
                'haplotypes': args.haplotypes or [],
                'save_fig': args.output,
                'width': args.width,
                'height': args.height,
            }
        }
    elif plot_type == 'gene':
        config = {
            'data': {
                'gff_file': args.gff,
                'genes': args.genes,
                'toml_file': args.toml
            },
            'plot': {
                'upstream': args.upstream,
                'downstream': args.downstream,
                'width': args.width,
                'height': args.height,
                'save_fig': args.output
            }
        }
    elif plot_type == 'network':
        config = {
            'data': {
                'file': args.file,
                'sample_hap': args.sample_hap,
                'sample_group': args.sample_group,
                'group_index': args.group_index
            },
            'plot': {
                'width': args.width,
                'height': args.height,
                'save_fig': args.output,
                'show_node_label': getattr(args, 'show_node_label', False),
                'node_label_size': getattr(args, 'node_label_size', 10),
                'node_label_color': getattr(args, 'node_label_color', 'black'),
                'show_edge_label': getattr(args, 'show_edge_label', True),
                'edge_label_symbol': getattr(args, 'edge_label_symbol', None),
                'edge_label_size': getattr(args, 'edge_label_size', 8),
                'edge_label_color': getattr(args, 'edge_label_color', 'black'),
                'layout': getattr(args, 'layout', 'spring'),
                'k': getattr(args, 'k', None),
                'node_scale_factor': getattr(args, 'node_scale_factor', 1.0),
                'node_color': getattr(args, 'node_color', None),
                'node_alpha': getattr(args, 'node_alpha', 0.8),
                'edge_color': getattr(args, 'edge_color', 'gray'),
                'edge_width': getattr(args, 'edge_width', 1.0),
                'edge_alpha': getattr(args, 'edge_alpha', 0.6),
                'edge_style': getattr(args, 'edge_style', 'line'),
                'seed': getattr(args, 'seed', None),
                'title': getattr(args, 'title', 'Haplotype Network')
            }
        }
    else:
        logger.error(f'Unknown plot type: {plot_type}')
        return

    try:
        plot_instance = manager.use(plot_type, config)
        plot_instance.plot()
        logger.info(f'Successfully created {plot_type} plot: {config["plot"]["save_fig"]}')
    except ValueError as e:
        logger.error(f'Error creating {plot_type} plot: {e}')
    except Exception as e:
        logger.error(f'Unexpected error creating {plot_type} plot: {e}')

    end = time.time()
    logger.info("plot function runs {:.2f} seconds".format(end - start))
