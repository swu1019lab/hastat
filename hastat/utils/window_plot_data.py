"""Validate existing window results for plotting; never recalculate statistics."""
import csv
import numpy as np
import pandas as pd
from hastat.utils.sample_data import numeric_trait


def read_window_track(path, statistic):
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
    frame[statistic] = numeric_trait(frame[statistic])
    if statistic == 'pi' and (frame.pi < 0).any():
        raise ValueError('Nucleotide diversity pi cannot be negative')
    keys = [name for name in ('chrom', 'groups', 'start', 'end') if name in frame]
    if frame.duplicated(keys).any():
        raise ValueError('Duplicate windows within a chromosome/group are ambiguous')
    return frame
