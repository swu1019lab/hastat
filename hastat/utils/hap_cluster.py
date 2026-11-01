"""Hierarchical display clustering of discrete haplotype profiles."""
import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial import distance

from hastat.log.logger import logger


def cluster_haplotypes(table, haplotypes, method='average'):
    """Return ordered IDs and linkage on common complete loci (Hamming distance).

    Stable ID order resolves input-order dependence in distance ties. Each profile
    has equal weight; sample counts affect labels only. This is a display tree,
    not a phylogenetic inference or statistical significance test.
    """
    allowed = ('single', 'complete', 'average', 'weighted')
    if method not in allowed:
        raise ValueError('hap_cluster_linkage must be single, complete, average or weighted '
                         '(Hamming distance); ward, centroid and median require Euclidean distances')
    names = sorted(haplotypes)
    if len(set(names)) != len(names):
        raise ValueError('Clustering requires unique haplotype IDs')
    values = table[names].to_numpy(dtype=float)
    if not np.isin(values, [-1, 0, 1, 2]).all():
        raise ValueError('Clustering genotypes must be -1 (missing), 0, 1 or 2')
    if len(names) < 2:
        return names, None
    complete = (values != -1).all(axis=1)
    excluded = int((~complete).sum())
    if excluded:
        logger.warning(f'Haplotype clustering excludes {excluded} loci with missing calls; '
                       'all pairs use the same complete loci')
    if not complete.any():
        raise ValueError('No complete loci shared by displayed haplotypes for clustering')
    distances = distance.pdist(values[complete].T, metric='hamming')
    linkage = hierarchy.linkage(distances, method=method)
    leaves = hierarchy.leaves_list(linkage)
    return [names[i] for i in leaves], linkage
