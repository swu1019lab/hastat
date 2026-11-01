"""Shared interval, genotype-profile and sample-count operations."""
from collections import defaultdict

import numpy as np
import pandas as pd

from hastat.log.logger import logger


def normalize_loci(loci):
    loci = list(loci)
    if not loci:
        raise ValueError('At least one genomic interval is required')
    by_chrom = defaultdict(list)
    for chrom, start, end in loci:
        if not isinstance(chrom, str) or not chrom:
            raise ValueError('Chromosome names must be non-empty strings')
        if (not isinstance(start, (int, np.integer)) or not isinstance(end, (int, np.integer))
                or start < 1 or end < start):
            raise ValueError('Regions must be 1-based inclusive with 1 <= start <= end')
        by_chrom[chrom].append((start, end))
    merged = []
    for chrom in sorted(by_chrom):
        intervals = []
        for start, end in sorted(by_chrom[chrom]):
            if intervals and start <= intervals[-1][1] + 1:
                intervals[-1] = (intervals[-1][0], max(end, intervals[-1][1]))
            else:
                intervals.append((start, end))
        merged.extend((chrom, start, end) for start, end in intervals)
    return merged


def haplotype_outputs(data, prefix='Hap'):
    """Build both outputs from the same complete, canonically ordered profiles."""
    if data.index.has_duplicates or not data.columns.is_unique:
        raise ValueError('Haplotype input must have unique variant records and sample names')
    data = data.sort_index()
    if len(data) == 0:
        return (pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt']),
                pd.DataFrame(columns=['samples', 'haplotypes']))
    if data.shape[1] == 0 or not data.isin([0, 1, 2]).all().all():
        raise ValueError('Haplotype grouping requires samples with complete 0/1/2 states')
    hap, inverse, counts = np.unique(data.T.to_numpy(dtype=np.int8), axis=0,
                                     return_inverse=True, return_counts=True)
    order = sorted(range(len(hap)), key=lambda i: (int(counts[i]), tuple(hap[i])), reverse=True)
    ranks = np.empty(len(order), dtype=int)
    ranks[order] = np.arange(len(order))
    labels = [f'{prefix}{i + 1}' for i in range(len(order))]
    table = pd.DataFrame(hap[order].T, index=data.index, columns=labels).reset_index()
    groups = pd.DataFrame({'samples': data.columns.tolist(),
                           'haplotypes': [labels[ranks[i]] for i in inverse]})
    return table, groups


def count_haplotypes(assignments, sample_groups, available_samples):
    """Count assigned individuals, optionally restricted to a validated group file."""
    order = sorted(assignments.haplotypes.unique(), key=lambda hap: int(hap[3:]))
    if sample_groups is None:
        counts = assignments.groupby('haplotypes').size().reindex(order, fill_value=0)
        counts.index.name = 'haplotypes'
        return counts.astype('int64').rename('samples').to_frame()
    if (not isinstance(sample_groups, pd.DataFrame) or len(sample_groups.columns) != 2
            or set(sample_groups.columns) != {'samples', 'groups'}):
        raise ValueError('Frequency groups must contain exactly samples and groups columns')
    metadata = sample_groups.loc[:, ['samples', 'groups']].copy()
    if metadata.empty or any(not isinstance(value, str) or not value.strip()
                             for value in metadata.to_numpy().flat):
        raise ValueError('Frequency groups require non-empty string sample and group names')
    if metadata.samples.duplicated().any():
        raise ValueError('Frequency group file must assign each sample exactly once')
    unknown = metadata.loc[~metadata.samples.isin(available_samples), 'samples'].tolist()
    if unknown:
        raise ValueError(f"Samples not found in VCF: {', '.join(unknown)}")
    labels = sorted(metadata.groups.unique())
    if 'haplotypes' in labels:
        raise ValueError("Population label 'haplotypes' is reserved for the frequency row identifier")
    if len(metadata) < len(available_samples):
        logger.info(f'Frequency counts use the {len(metadata)} listed samples; '
                    'haplotype definitions and ranks still use all VCF samples')
    joined = assignments.merge(metadata, on='samples', how='inner', validate='one_to_one')
    if joined.empty:
        return pd.DataFrame(index=pd.Index([], name='haplotypes'), columns=labels, dtype='int64')
    counts = joined.groupby(['haplotypes', 'groups']).size().unstack(fill_value=0)
    observed = [hap for hap in order if hap in counts.index]
    counts = counts.reindex(index=observed, columns=labels, fill_value=0).astype('int64')
    counts.index.name, counts.columns.name = 'haplotypes', None
    return counts
