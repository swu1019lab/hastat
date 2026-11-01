"""Shared sample CSV parsing; identifiers and labels are literal text."""
import csv
import numpy as np
import pandas as pd


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
    frame = sample_frame(frame, 'haplotypes' if label == 'haplotypes' else 'groups')
    frame = frame.rename(columns={'sample': 'samples'})
    if frame.empty or frame[label].nunique() < min_groups:
        raise ValueError(f'Assignment file requires at least {min_groups} non-empty group(s)')
    if available_samples is not None:
        unknown = frame.loc[~frame.samples.isin(available_samples), 'samples'].tolist()
        if unknown:
            raise ValueError(f"Samples not found in VCF: {', '.join(unknown)}")
    return frame


def sample_frame(frame, kind):
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


def numeric_trait(values):
    # Missing tokens apply only to measurements, never to sample/group labels.
    missing = {'', 'NA', 'N/A', 'NAN', 'NULL', 'NONE', '<NA>', '#N/A', '#NA',
               '#N/A N/A', '-NAN', '1.#QNAN', '-1.#IND', '1.#IND', '-1.#QNAN'}
    values = values.map(lambda x: np.nan if isinstance(x, str) and x.strip().upper() in missing else x)
    result = pd.to_numeric(values, errors='raise').astype(float)
    return result.where(np.isfinite(result), np.nan)


