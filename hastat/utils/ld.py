"""Rogers-Huff r-squared matrices; see LD_METHODS.md for scope and missingness."""
import allel
import numpy as np
from scipy.spatial.distance import squareform


def ld_r_squared(genotypes):
    """Return n_variants x n_variants r², preserving undefined values as NaN.

    Rows are variants, columns are diploid samples. Values are 0/1/2 alternate
    allele counts, with -1 or NaN for missing genotypes. No imputation is used.
    The library estimates off-diagonal r from pairwise complete observations.
    """
    data = np.asarray(genotypes)
    if data.ndim != 2:
        raise ValueError('LD genotypes must be a two-dimensional variants-by-samples array')
    try:
        data = data.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError('LD genotypes must be 0/1/2, -1 or NaN') from exc
    if not np.all(np.isnan(data) | np.isin(data, [-1, 0, 1, 2])):
        raise ValueError('LD genotypes must be 0/1/2, -1 or NaN')
    data = np.ascontiguousarray(np.where(np.isnan(data), -1, data), dtype=np.int8)
    n_variants, n_samples = data.shape
    if n_variants == 0 or n_samples < 2:
        return np.full((n_variants, n_variants), np.nan, dtype=np.float32)
    if n_variants > 1:
        # scikit-allel returns a scalar for exactly two variants.
        r = np.atleast_1d(allel.rogers_huff_r(data))
        if r.shape != (n_variants * (n_variants - 1) // 2,):
            raise RuntimeError('Unexpected condensed LD shape from scikit-allel')
        squared = np.square(r)
        squared[~np.isfinite(squared)] = np.nan
        # float32 arithmetic can give r² slightly above one; NaN stays NaN.
        result = squareform(np.clip(squared, 0, 1))
    else:
        result = np.full((1, 1), np.nan, dtype=np.float32)
    called = data >= 0
    valid = ((called.sum(axis=1) >= 2)
             & (np.where(called, data, 3).min(axis=1)
                < np.where(called, data, -1).max(axis=1)))
    # Allele polymorphism alone is insufficient: all 0/1 calls have zero dosage variance.
    np.fill_diagonal(result, np.where(valid, 1.0, np.nan))
    return result
