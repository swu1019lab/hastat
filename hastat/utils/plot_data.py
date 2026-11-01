"""Validated sample joins for box, bar and pie plots; no output sidecars."""
from hastat.log.logger import logger
from hastat.utils.sample_data import read_sample_csv, sample_frame, numeric_trait


def select_haplotypes(frame, labels):
    if labels is None or len(labels) == 0:
        return frame.copy()
    labels = [str(x).strip() for x in labels]
    if set(labels) - set(frame['haplotypes']):
        raise ValueError('Requested haplotype is absent from the haplotype table')
    return frame.loc[frame['haplotypes'].isin(labels)].copy()


def merge_samples(haplotypes, metadata):
    merged = haplotypes.merge(metadata, on='sample', how='inner', validate='one_to_one')
    if merged.empty:
        raise ValueError('No overlapping samples to plot')
    left = len(haplotypes) - len(merged)
    right = len(metadata) - len(merged)
    if left or right:
        logger.warning(f'Sample join excluded {left} haplotype rows and {right} metadata rows without a match')
    return merged


def grouped_counts(data, name, calc_percentage=False):
    if name in ('sample', 'haplotypes') or name not in data:
        raise ValueError('Select a metadata group column, not the sample or haplotype column')
    counts = data.groupby(['haplotypes', name], observed=True).size().unstack(fill_value=0)
    if counts.empty:
        raise ValueError('No grouped samples to plot')
    return counts.div(counts.sum()) if calc_percentage else counts
