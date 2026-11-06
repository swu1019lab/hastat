"""Read, validate and prepare sample and window tables."""
import csv
import numpy as np
import pandas as pd
from hastat.log.logger import logger


def read_sample_csv(path):
    # Inspect the actual header before pandas can rename duplicate columns.
    with open(path, encoding='utf-8-sig', newline='') as handle:
        header = next(csv.reader(handle), [])
    if len(header) < 2 or len(set(header)) != len(header) or any(not x.strip() for x in header):
        raise ValueError('Sample CSV requires at least two uniquely named, nonempty columns')
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def validate_sample_ids(frame, column='samples'):
    if not isinstance(frame, pd.DataFrame) or not frame.columns.is_unique or column not in frame:
        raise ValueError(f'Sample table requires unique columns including {column}')
    result = frame.copy()
    if result[column].isna().any():
        raise ValueError('Missing sample identifiers')
    result[column] = result[column].astype(str).str.strip()
    if result[column].eq('').any() or result[column].duplicated().any():
        raise ValueError('Sample table requires one row per unique sample; empty or duplicate ID found')
    return result


def read_assignments(path, label='groups', available_samples=None, min_groups=1):
    """Read a two-column assignment file with unique samples and nonempty labels."""
    frame = read_sample_csv(path)
    if len(frame.columns) != 2:
        raise ValueError('Assignment file requires exactly two columns: sample and group')
    frame.columns = ['samples', label]
    frame = normalize_sample_table(frame, 'haplotypes' if label == 'haplotypes' else 'groups')
    frame = frame.rename(columns={'sample': 'samples'})
    if frame.empty or frame[label].nunique() < min_groups:
        raise ValueError(f'Assignment file requires at least {min_groups} non-empty group(s)')
    if available_samples is not None:
        unknown = frame.loc[~frame.samples.isin(available_samples), 'samples'].tolist()
        if unknown:
            raise ValueError(f"Samples not found in VCF: {', '.join(unknown)}")
    return frame


def normalize_sample_table(frame, kind):
    if not isinstance(frame, pd.DataFrame):
        raise ValueError('Input must be a DataFrame')
    if len(frame.columns) < 2 or not frame.columns.is_unique:
        raise ValueError('Sample table requires at least two unique columns')
    if kind == 'haplotypes' and len(frame.columns) != 2:
        raise ValueError('Haplotype table requires exactly two columns')
    result = frame.copy()
    result.columns = (['sample', 'haplotypes'] if kind == 'haplotypes'
                      else ['sample'] + list(result.columns[1:]))
    if not result.columns.is_unique or (kind != 'haplotypes' and 'haplotypes' in result):
        raise ValueError('Metadata columns cannot be named sample or haplotypes')
    result = validate_sample_ids(result, 'sample')
    text_columns = result.columns[1:] if kind in ('haplotypes', 'groups') else []
    for column in text_columns:
        if result[column].isna().any():
            raise ValueError(f'Missing values in {column}')
        result[column] = result[column].astype(str).str.strip()
        if result[column].eq('').any():
            raise ValueError(f'Empty values in {column}')
    return result


def parse_numeric_values(values):
    # Missing tokens apply only to measurements, never to sample/group labels.
    missing = {'', 'NA', 'N/A', 'NAN', 'NULL', 'NONE', '<NA>', '#N/A', '#NA',
               '#N/A N/A', '-NAN', '1.#QNAN', '-1.#IND', '1.#IND', '-1.#QNAN'}
    values = values.map(lambda x: np.nan if isinstance(x, str) and x.strip().upper() in missing else x)
    result = pd.to_numeric(values, errors='raise').astype(float)
    return result.where(np.isfinite(result), np.nan)


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


def read_window_table(path, statistic):
    if statistic not in ('pi', 'fst', 'tajima_d'):
        raise ValueError('Window track statistic must be pi, fst or tajima_d')
    with open(path, encoding='utf-8-sig', newline='') as handle:
        header = next(csv.reader(handle), [])
    if len(header) != len(set(header)) or any(not name.strip() for name in header):
        raise ValueError('Window CSV requires uniquely named, nonempty columns')
    required = {'start', 'end', statistic}
    if not required.issubset(header):
        raise ValueError(f'Window CSV requires columns: start, end, {statistic}')
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    for label in ('chrom', 'groups'):
        if label in frame:
            frame[label] = frame[label].str.strip()
            if frame[label].eq('').any():
                raise ValueError(f'Window CSV contains empty {label} labels')
    for coord in ('start', 'end'):
        values = pd.to_numeric(frame[coord], errors='raise')
        if (not np.isfinite(values).all() or (values < 1).any()
                or (values != np.floor(values)).any()):
            raise ValueError('Window coordinates must be positive integers (1-based inclusive)')
        frame[coord] = values.astype('int64')
    if (frame.end < frame.start).any():
        raise ValueError('Window end must be greater than or equal to start')
    raw = frame[statistic].str.strip().str.lower()
    if raw.str.fullmatch(r'[+-]?inf(?:inity)?').any():
        raise ValueError('Window statistics cannot contain infinity; use NaN for undefined values')
    frame[statistic] = parse_numeric_values(frame[statistic])
    if statistic == 'pi' and (frame.pi < 0).any():
        raise ValueError('Nucleotide diversity pi cannot be negative')
    keys = [name for name in ('chrom', 'groups', 'start', 'end') if name in frame]
    if frame.duplicated(keys).any():
        raise ValueError('Duplicate windows within a chromosome/group are ambiguous')
    return frame
