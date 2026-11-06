"""Population genetic statistics on diploid variant data."""
import allel
import numpy as np
from scipy.spatial.distance import squareform
from hastat.log.logger import logger
from hastat.utils.variant import VariantData, normalize_loci

TAJIMA_COLUMNS = ['chrom', 'start', 'end', 'counts', 'tajima_d']


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


def _validate_windows(size, step):
    if any(not isinstance(x, (int, np.integer)) or x < 1 for x in (size, step)):
        raise ValueError("Window size and step must be positive integers")


def _population_windows(start, end, size, step):
    VariantData._validate_region(start, end)
    _validate_windows(size, step)
    windows = []
    for left in range(start, end + 1, step):
        right = min(left + size - 1, end)
        windows.append((left, right))
        if right == end:
            break
    return np.asarray(windows, dtype=int)


def _population_allele_counts(variants, chrom, start, end, sample_list, het_threshold):
    """Count chromosomes, not dosage categories, at biallelic SNPs."""
    data = variants.extract_genotypes(chrom, start, end, 'code2', het_threshold=het_threshold)
    samples = list(data.columns) if sample_list is None else list(sample_list)
    if not samples or len(set(samples)) != len(samples):
        raise ValueError("Population sample list must be non-empty and unique")
    unknown = [s for s in samples if s not in data.columns]
    if unknown:
        raise ValueError(f"Samples not found in VCF: {', '.join(map(str, unknown))}")
    data = data.loc[:, samples]
    is_snp = [ref in ('A', 'C', 'G', 'T') and alt in ('A', 'C', 'G', 'T')
              and ref != alt for _, _, ref, alt in data.index]
    if not all(is_snp):
        logger.warning("Pi/Tajima's D use biallelic SNPs only; excluding non-SNP records")
    data = data.loc[is_snp]
    pos = data.index.get_level_values('pos').to_numpy(dtype=int)
    if len(np.unique(pos)) != len(pos):
        raise ValueError("Pi/Tajima's D require one biallelic SNP record per position")
    if het_threshold is not None:
        logger.warning("Heterozygosity filtering changes the site frequency spectrum; "
                       "these estimates describe retained SNPs only")
    calls = np.asarray(data.to_numpy().tolist(), dtype=np.int8).reshape(len(data), len(samples), 2)
    ac = allel.GenotypeArray(calls).count_alleles(max_allele=1)
    return pos, ac, 2 * len(samples)


def calculate_pi(variants, chrom_name, chrom_start, chrom_end, size=1, step=1, sample_list=None, het_threshold=None):
    """
    Get pi data of target samples using sliding window method

    :param chrom_name: chromosome name
    :param chrom_start: chromosome start position
    :param chrom_end: chromosome end position
    :param size: window size
    :param step: step size
    :param sample_list: a list contained sample names
    :param het_threshold: heterozygosity rate threshold for filtering variants
    :return: a numpy array contained pi data
    """
    windows = _population_windows(chrom_start, chrom_end, size, step)
    pos, ac, _ = _population_allele_counts(variants,
        chrom_name, chrom_start, chrom_end, sample_list, het_threshold)
    pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac, windows=windows)
    # scikit-allel fills sites with fewer than two called alleles with zero.
    # Such a recorded SNP is unknown, not known monomorphic: invalidate
    # its windows instead of silently reducing diversity.
    insufficient, _, _ = allel.windowed_statistic(
        pos, ac.sum(axis=1) < 2, statistic=np.any, windows=windows, fill=False)
    if np.any(insufficient):
        pi[insufficient] = np.nan
        logger.warning(f"Pi is NaN in {np.count_nonzero(insufficient)} windows with uncalled SNPs")
    return np.column_stack([[chrom_name] * len(windows), windows, n_bases, counts, pi])


def calculate_multi_loci_pi(variants, loci, size=1, step=1, sample_list=None, het_threshold=None):
    """
    Get pi data from multiple regions (for homologous genes)

    :param loci: list of tuples containing (chrom_name, chrom_start, chrom_end)
    :param size: window size
    :param step: step size
    :param sample_list: a list contained sample names
    :param het_threshold: heterozygosity rate threshold for filtering variants
    :return: a numpy array contained pi data
    """
    logger.info(f"Calculating Pi data from {len(loci)} regions")

    all_pi_data = []
    for chrom_name, chrom_start, chrom_end in loci:
        pi_data = calculate_pi(variants, chrom_name, chrom_start, chrom_end, size, step, sample_list, het_threshold)
        if pi_data.size > 0:
            all_pi_data.append(pi_data)

    if not all_pi_data:
        return np.array([])

    return np.vstack(all_pi_data)


def calculate_fst(variants, chrom_name, chrom_start, chrom_end, sample_list1: list, sample_list2: list, size=1,
                 step=1, het_threshold=None):
    """
    Get Fst data of target samples using sliding window method

    :param chrom_name: chromosome name
    :param chrom_start: chromosome start position
    :param chrom_end: chromosome end position
    :param sample_list1: a list contained sample names
    :param sample_list2: a list contained sample names
    :param size: sliding window size
    :param step: step size
    :param het_threshold: heterozygosity rate threshold for filtering variants
    :return: a numpy array contained Fst data
    """
    windows = _population_windows(chrom_start, chrom_end, size, step)
    groups = []
    available = set(variants.vcf_in.header.samples)
    for samples in (sample_list1, sample_list2):
        if samples is None or isinstance(samples, (str, bytes)):
            raise ValueError("Fst populations must be non-empty sample lists")
        samples = list(samples)
        if not samples or len(set(samples)) != len(samples):
            raise ValueError("Fst population sample lists must be non-empty and unique")
        unknown = [s for s in samples if s not in available]
        if unknown:
            raise ValueError(f"Samples not found in VCF: {', '.join(map(str, unknown))}")
        groups.append(samples)
    if set(groups[0]) & set(groups[1]):
        raise ValueError("Fst populations must not share samples")
    logger.info(f"Calculating Fst between {len(groups[0])} and {len(groups[1])} samples")
    geno_data = variants.extract_genotypes(chrom_name, chrom_start, chrom_end, 'code2', het_threshold=het_threshold)
    selected = groups[0] + groups[1]
    split = len(groups[0])
    sub_pops = [list(range(split)), list(range(split, len(selected)))]
    pos = geno_data.index.get_level_values('pos').to_numpy(dtype=int)
    # Preserve counts as input biallelic-record counts, including records
    # without enough called individuals for variance-component estimation.
    counts = (np.searchsorted(pos, windows[:, 1], side='right') -
              np.searchsorted(pos, windows[:, 0], side='left'))
    fst = np.full(len(windows), np.nan)
    if len(geno_data):
        calls = np.asarray(geno_data.loc[:, selected].to_numpy().tolist(), dtype=np.int8)
        called = np.all(calls >= 0, axis=2)
        n1, n2 = called[:, :split].sum(axis=1), called[:, split:].sum(axis=1)
        # For two populations, WC84 requires both represented, n_bar > 1
        # and n_C > 0. One called individual in a group is not by itself
        # invalid when the other group has at least two.
        usable = (n1 > 0) & (n2 > 0) & (n1 + n2 > 2)
        excluded = int(np.count_nonzero(~usable))
        if excluded:
            logger.warning(f"Fst excludes {excluded} records with insufficient called samples "
                           "for Weir-Cockerham variance components")
        if np.any(usable):
            with np.errstate(divide='ignore', invalid='ignore'):
                fst, _, _ = allel.windowed_weir_cockerham_fst(
                    pos[usable], allel.GenotypeArray(calls[usable]), sub_pops,
                    windows=windows, max_allele=1)
            fst[~np.isfinite(fst)] = np.nan
    return np.column_stack([[chrom_name] * len(windows), windows, counts, fst])


def calculate_multi_loci_fst(variants, loci, sample_list1: list, sample_list2: list, size=1, step=1, het_threshold=None):
    """
    Get Fst data from multiple regions (for homologous genes)

    :param loci: list of tuples containing (chrom_name, chrom_start, chrom_end)
    :param sample_list1: a list contained sample names
    :param sample_list2: a list contained sample names
    :param size: sliding window size
    :param step: step size
    :param het_threshold: heterozygosity rate threshold for filtering variants
    :return: a numpy array contained Fst data
    """
    logger.info(f"Calculating Fst data from {len(loci)} regions")

    all_fst_data = []
    for chrom_name, chrom_start, chrom_end in loci:
        fst_data = calculate_fst(variants, chrom_name, chrom_start, chrom_end, sample_list1, sample_list2, size, step, het_threshold)
        if fst_data.size > 0:
            all_fst_data.append(fst_data)

    if not all_fst_data:
        return np.array([])

    return np.vstack(all_fst_data)


def calculate_tajima_d(variants, chrom_name, chrom_start, chrom_end, size=1, step=1, sample_list=None, het_threshold=None):
    """
    Get Tajima's D data of target samples using sliding window method

    :param chrom_name: chromosome name
    :param chrom_start: chromosome start position
    :param chrom_end: chromosome end position
    :param size: window size
    :param step: step size
    :param sample_list: a list contained sample names
    :param het_threshold: heterozygosity rate threshold for filtering variants
    :return: a numpy array with TAJIMA_COLUMNS

    Keep the selected population fixed and exclude each SNP with any
    missing diploid call in that  Apply the original Tajima D
    formula to the remaining complete sites; no variable-n approximation.
    """
    windows = _population_windows(chrom_start, chrom_end, size, step)
    pos, ac, n_chromosomes = _population_allele_counts(variants,
        chrom_name, chrom_start, chrom_end, sample_list, het_threshold)

    complete = ac.sum(axis=1) == n_chromosomes
    rows = []
    for left, right in windows:
        first, last = np.searchsorted(pos, left), np.searchsorted(pos, right, side='right')
        window_ac = ac[first:last][complete[first:last]]
        counts = int(last - first)
        segregating = int(window_ac.count_segregating())
        value = np.nan
        if n_chromosomes >= 4 and segregating >= 3:
            value = float(allel.tajima_d(window_ac, min_sites=3))
            if not np.isfinite(value):
                value = np.nan
        rows.append([chrom_name, int(left), int(right), counts, value])
    if np.any(~complete):
        logger.warning(f"Tajima's D excludes {np.count_nonzero(~complete)} SNPs with missing "
                       "calls in the selected population; uses complete sites with a fixed "
                       "sample set. Non-random missingness can bias the retained spectrum.")
    return np.asarray(rows, dtype=str)


def calculate_multi_loci_tajima_d(variants, loci, size=1, step=1, sample_list=None, het_threshold=None):
    """
    Get Tajima's D data from multiple regions (for homologous genes)

    :param loci: list of tuples containing (chrom_name, chrom_start, chrom_end)
    :param size: window size
    :param step: step size
    :param sample_list: a list contained sample names
    :param het_threshold: heterozygosity rate threshold for filtering variants
    :return: a numpy array contained Tajima's D data
    """
    logger.info(f"Calculating Tajima's D data from {len(loci)} regions")

    all_tajima_d_data = []
    for chrom_name, chrom_start, chrom_end in loci:
        tajima_d_data = calculate_tajima_d(variants, chrom_name, chrom_start, chrom_end, size, step, sample_list, het_threshold)
        if tajima_d_data.size > 0:
            all_tajima_d_data.append(tajima_d_data)

    if not all_tajima_d_data:
        return np.empty((0, len(TAJIMA_COLUMNS)), dtype=str)

    return np.vstack(all_tajima_d_data)


def calculate_ld(variants, chrom_name, chrom_start, chrom_end):
    """Return a square r² matrix in VCF record order (empty regions: 0 x 0)."""
    logger.info(f"Calculating LD data in target region {chrom_name}:{chrom_start}-{chrom_end}")
    geno_data = variants.extract_genotypes(chrom_name, chrom_start, chrom_end, 'code1')
    if not geno_data.index.is_unique:
        raise ValueError('Duplicate variant records in LD input')
    return ld_r_squared(geno_data.to_numpy())


def calculate_multi_loci_ld(variants, loci):
    """Return r² for sorted, merged intervals; overlapping requests count once.

    Matrix order follows extract_multi_loci_genotypes(normalize_loci(loci)).
    Distinct chromosomes are allowed, but r² is not evidence of physical linkage.
    """
    loci = normalize_loci(loci)
    logger.info(f"Calculating LD data from {len(loci)} regions")
    geno_data = variants.extract_multi_loci_genotypes(loci, 'code1')
    if not geno_data.index.is_unique:
        raise ValueError('Duplicate variant records in LD input')
    return ld_r_squared(geno_data.to_numpy())


def calculate_heterozygosity(variants, chrom_name, chrom_start, chrom_end):
    """
    Get heterozygosity data of all samples

    :return: a numpy array contained heterozygosity data
    """
    geno_data = variants.extract_genotypes(chrom_name, chrom_start, chrom_end, 'code2')
    g = allel.GenotypeArray(geno_data.to_numpy().tolist())
    # the rate of observed heterozygosity
    het_o = allel.heterozygosity_observed(g)
    # the expected rate of heterozygosity
    af = g.count_alleles().to_frequencies()
    het_e = allel.heterozygosity_expected(af, ploidy=2)
    return np.column_stack([het_o, het_e])


def calculate_inbreeding_coefficient(variants, chrom_name, chrom_start, chrom_end):
    """
    Get inbreeding coefficient of all samples

    :return: a numpy array contained inbreeding coefficient
    """
    geno_data = variants.extract_genotypes(chrom_name, chrom_start, chrom_end, 'code2')
    g = allel.GenotypeArray(geno_data.to_numpy().tolist())
    return allel.inbreeding_coefficient(g)
