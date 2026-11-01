# -*- coding: utf-8 -*-
# @Time    : 2024/1/12 17:32
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : gene.py

import numpy as np
import pandas as pd
import networkx as nx
from collections import defaultdict
import allel
from hastat.log.logger import logger
from hastat.utils.haplotype import normalize_loci, haplotype_outputs, count_haplotypes
from hastat.utils.ld import ld_r_squared


TAJIMA_COLUMNS = ['chrom', 'start', 'end', 'counts', 'tajima_d']
TAJIMA_QC_COLUMNS = TAJIMA_COLUMNS + [
    'n_samples', 'n_chromosomes', 'sites_used', 'sites_excluded_missing',
    'segregating_sites',
]


class GeneFeature(object):
    """
    A class for processing gffutils.FeatureDB object
    """

    def __init__(self, gff_in) -> None:
        """
        :param gff_in: gffutils.FeatureDB object
        """
        self.gff_in = gff_in

    def get_locus(self, gene_id, upstream=0, downstream=0):
        """
        Get the locus of the gene

        :param gene_id: gene id
        :param upstream: the length of upstream with respect to the gene, default 0 bp
        :param downstream: the length of downstream with respect to the gene, default 0 bp
        :return: a tuple containing the locus of the gene
        """
        gene = self.gff_in[gene_id]
        if upstream < 0 or downstream < 0:
            raise ValueError("Upstream and downstream lengths must be non-negative")
        if gene.strand == "+":
            start = max(1, gene.start - upstream)
            end = gene.end + downstream
        else:
            start = max(1, gene.start - downstream)
            end = gene.end + upstream
        return gene.seqid, start, end

    def get_merged_locus(self, gene_ids, upstream=0, downstream=0):
        """
        Get the loci of multiple homologous genes (supports cross-chromosome genes)

        :param gene_ids: list of gene ids
        :param upstream: the length of upstream with respect to the genes, default 0 bp
        :param downstream: the length of downstream with respect to the genes, default 0 bp
        :return: a list of tuples containing the loci of the genes [(chrom1, start1, end1), ...]
        """
        if not gene_ids:
            raise ValueError("Gene IDs list cannot be empty")
        
        logger.info(f"Processing {len(gene_ids)} homologous genes: {', '.join(gene_ids)}")
        
        gene_loci = []
        chromosomes = set()
        
        # Get locus for each gene
        for gene_id in gene_ids:
            try:
                chrom, start, end = self.get_locus(gene_id, upstream, downstream)
                gene_loci.append((chrom, start, end))
                chromosomes.add(chrom)
                logger.info(f"Gene {gene_id}: {chrom}:{start}-{end}")
            except KeyError:
                logger.error(f"Gene ID '{gene_id}' not found in GFF file")
                raise ValueError(f"Gene ID '{gene_id}' not found in GFF file")
        
        # Group genes by chromosome and merge overlapping regions
        from collections import defaultdict
        chrom_regions = defaultdict(list)
        
        for chrom, start, end in gene_loci:
            chrom_regions[chrom].append((start, end))
        
        # Merge overlapping regions within each chromosome
        merged_loci = []
        for chrom, regions in chrom_regions.items():
            # Sort regions by start position
            regions.sort()
            merged_regions = []
            
            for start, end in regions:
                if merged_regions and start <= merged_regions[-1][1]:
                    # Overlapping or adjacent regions, merge them
                    merged_regions[-1] = (merged_regions[-1][0], max(merged_regions[-1][1], end))
                else:
                    # Non-overlapping region, add as new
                    merged_regions.append((start, end))
            
            # Add merged regions for this chromosome
            for start, end in merged_regions:
                merged_loci.append((chrom, start, end))
        
        logger.info(f"Found genes on {len(chromosomes)} chromosome(s): {', '.join(sorted(chromosomes))}")
        for chrom, start, end in merged_loci:
            logger.info(f"Merged region: {chrom}:{start}-{end}")
        
        return merged_loci

    def get_feature_data(self, gene_id, feature_type='CDS'):
        logger.info("Retrieving {} data of {} from FeatureDB".format(feature_type, gene_id))
        gene = self.gff_in[gene_id]
        return np.array([[f.start, f.end] for f in self.gff_in.region(
            seqid=gene.seqid,
            start=gene.start,
            end=gene.end,
            featuretype=feature_type
        )])


class GeneVariant(object):
    """
    A class for processing variants data of a gene or region
    """

    def __init__(self, vcf_in) -> None:
        """
        :param vcf_in: VariantFile object
        """
        self.vcf_in = vcf_in

    def _haplotype_result(self, loci, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
        loci = normalize_loci(loci)
        options = dict(het_threshold=het_threshold, maf_threshold=maf_threshold,
                       max_missing_threshold=max_missing_threshold)
        data = (self.get_geno_data(*loci[0], **options) if len(loci) == 1
                else self.get_multi_loci_geno_data(loci, **options))
        if data.index.has_duplicates:
            raise ValueError("Duplicate variant records in haplotype input")
        complete = data.isin([0, 1, 2]).all(axis=1)
        excluded = int((~complete).sum())
        if excluded:
            logger.warning(f"Haplotype grouping excludes {excluded} of {len(data)} records with missing calls; "
                           "all samples use the same complete loci")
        return haplotype_outputs(data.loc[complete], 'Hap')

    def hap_table(self, chrom_name, chrom_start, chrom_end, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
        """Return profiles on complete loci, retaining the original table columns."""
        return self._haplotype_result([(chrom_name, chrom_start, chrom_end)],
                                     het_threshold, maf_threshold, max_missing_threshold)[0]

    def multi_loci_hap_table(self, loci, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
        """Return profiles across the union of requested intervals."""
        return self._haplotype_result(loci, het_threshold, maf_threshold, max_missing_threshold)[0]

    def hap_groups(self, chrom_name, chrom_start, chrom_end, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
        """Assign samples using the exact same profiles/ranks as hap_table."""
        return self._haplotype_result([(chrom_name, chrom_start, chrom_end)],
                                     het_threshold, maf_threshold, max_missing_threshold)[1]

    def multi_loci_hap_groups(self, loci, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
        """Assign samples across the union of requested intervals."""
        return self._haplotype_result(loci, het_threshold, maf_threshold, max_missing_threshold)[1]

    def calc_hap_freq(self, chrom_name, chrom_start, chrom_end, sample_groups=None, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
        """Return integer sample counts per profile; group-file subsets affect counts only."""
        groups = self.hap_groups(chrom_name, chrom_start, chrom_end, het_threshold,
                                 maf_threshold, max_missing_threshold)
        return count_haplotypes(groups, sample_groups, list(self.vcf_in.header.samples))

    def calc_multi_loci_hap_freq(self, loci, sample_groups=None, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
        """Count multi-region profile assignments with the same rule as single-region freq."""
        groups = self.multi_loci_hap_groups(loci, het_threshold, maf_threshold, max_missing_threshold)
        return count_haplotypes(groups, sample_groups, list(self.vcf_in.header.samples))

    def get_geno_data(self, chrom_name, chrom_start, chrom_end, code='code1', het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
        """Read biallelic diploid calls whose POS lies in [start, end] (1-based).

        code1 is alternate-allele dosage (0/1/2); code2 contains two alleles.
        Any partially missing call becomes a wholly missing genotype (-1 or
        [-1, -1]). Filters use complete diploid calls across all VCF samples.
        """
        self._validate_region(chrom_start, chrom_end)
        if code not in ('code1', 'code2'):
            raise ValueError("Genotype code must be code1 or code2")
        for name, value, upper in [('maf_threshold', maf_threshold, 0.5),
                                   ('max_missing_threshold', max_missing_threshold, 1),
                                   ('het_threshold', het_threshold, 1)]:
            if value is not None and not 0 <= value <= upper:
                raise ValueError(f"{name} must be between 0 and {upper}")
        samples = list(self.vcf_in.header.samples)
        if not samples:
            raise ValueError("VCF must contain samples with diploid GT calls")
        indices, rows = [], []
        # pysam numeric coordinates are 0-based, half-open. GFF/VCF POS and
        # every public hastat region are 1-based, inclusive.
        for rec in self.vcf_in.fetch(chrom_name, chrom_start - 1, chrom_end):
            # fetch also returns deletions overlapping from before start.
            if not chrom_start <= rec.pos <= chrom_end:
                continue
            if rec.alleles is None or len(rec.alleles) != 2:
                continue
            if 'GT' not in rec.format:
                raise ValueError(f"Missing GT field at {rec.chrom}:{rec.pos}")
            pairs = []
            for sample in samples:
                gt = rec.samples[sample]['GT']
                if gt is None or gt == (None,):
                    gt = (None, None)
                if len(gt) != 2:
                    raise ValueError(f"Expected diploid GT at {rec.chrom}:{rec.pos}, sample {sample}")
                if any(a is not None and a not in (0, 1) for a in gt):
                    raise ValueError(f"Invalid biallelic GT at {rec.chrom}:{rec.pos}, sample {sample}")
                pairs.append([-1, -1] if None in gt else list(gt))
            dosage = [sum(gt) if gt[0] >= 0 else -1 for gt in pairs]
            called = [g for g in dosage if g >= 0]
            if (len(samples) - len(called)) / len(samples) > max_missing_threshold:
                continue
            if called:
                af = sum(called) / (2 * len(called))
                if min(af, 1 - af) < maf_threshold:
                    continue
            elif maf_threshold > 0:
                continue
            if het_threshold is not None:
                if not called or called.count(1) / len(called) > het_threshold:
                    continue
            indices.append((rec.chrom, rec.pos, rec.ref, rec.alts[0]))
            rows.append(dosage if code == 'code1' else pairs)
        index = pd.MultiIndex.from_tuples(indices, names=['chrom', 'pos', 'ref', 'alt'])
        columns = pd.Index(samples, name='samples')
        return pd.DataFrame(rows, index=index, columns=columns)

    @staticmethod
    def _validate_region(start, end):
        if (not isinstance(start, (int, np.integer)) or
                not isinstance(end, (int, np.integer)) or start < 1 or end < start):
            raise ValueError("Regions must be 1-based inclusive with 1 <= start <= end")

    @staticmethod
    def _validate_windows(size, step):
        if any(not isinstance(x, (int, np.integer)) or x < 1 for x in (size, step)):
            raise ValueError("Window size and step must be positive integers")

    @classmethod
    def _population_windows(cls, start, end, size, step):
        cls._validate_region(start, end)
        cls._validate_windows(size, step)
        windows = []
        for left in range(start, end + 1, step):
            right = min(left + size - 1, end)
            windows.append((left, right))
            if right == end:
                break
        return np.asarray(windows, dtype=int)

    def _population_allele_counts(self, chrom, start, end, sample_list, het_threshold):
        """Count chromosomes, not dosage categories, at biallelic SNPs."""
        data = self.get_geno_data(chrom, start, end, 'code2', het_threshold=het_threshold)
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

    def get_multi_loci_geno_data(self, loci, code='code1', het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
        """
        Get genotype data from multiple regions (for homologous genes)

        :param loci: list of tuples containing (chrom_name, chrom_start, chrom_end)
        :param code: genotype code, code1 or code2, default code1
        :param het_threshold: heterozygosity rate threshold for filtering variants
        :param maf_threshold: Minimum Minor Allele Frequency threshold
        :param max_missing_threshold: Maximum missing rate threshold
        :return: a DataFrame contained genotype data from all regions
        """
        logger.info(f"Extracting genotype data from {len(loci)} regions")
        
        all_data = []
        for chrom_name, chrom_start, chrom_end in loci:
            logger.info(f"Processing region: {chrom_name}:{chrom_start}-{chrom_end}")
            region_data = self.get_geno_data(chrom_name, chrom_start, chrom_end, code, het_threshold, maf_threshold, max_missing_threshold)
            if not region_data.empty:
                all_data.append(region_data)
        
        if not all_data:
            logger.warning("No genotype data found in any of the target regions")
            return pd.DataFrame(index=pd.MultiIndex.from_tuples([], names=['chrom', 'pos', 'ref', 'alt']),
                                columns=pd.Index(self.vcf_in.header.samples, name='samples'))
        
        # Concatenate all data
        combined_data = pd.concat(all_data, axis=0)
        logger.info(f"Combined data shape: {combined_data.shape}")
        
        return combined_data

    def get_pi_data(self, chrom_name, chrom_start, chrom_end, size=1, step=1, sample_list=None, het_threshold=None):
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
        windows = self._population_windows(chrom_start, chrom_end, size, step)
        pos, ac, _ = self._population_allele_counts(
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

    def get_multi_loci_pi_data(self, loci, size=1, step=1, sample_list=None, het_threshold=None):
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
            pi_data = self.get_pi_data(chrom_name, chrom_start, chrom_end, size, step, sample_list, het_threshold)
            if pi_data.size > 0:
                all_pi_data.append(pi_data)
        
        if not all_pi_data:
            return np.array([])
        
        return np.vstack(all_pi_data)

    def get_fst_data(self, chrom_name, chrom_start, chrom_end, sample_list1: list, sample_list2: list, size=1,
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
        windows = self._population_windows(chrom_start, chrom_end, size, step)
        groups = []
        available = set(self.vcf_in.header.samples)
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
        geno_data = self.get_geno_data(chrom_name, chrom_start, chrom_end, 'code2', het_threshold=het_threshold)
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

    def get_multi_loci_fst_data(self, loci, sample_list1: list, sample_list2: list, size=1, step=1, het_threshold=None):
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
            fst_data = self.get_fst_data(chrom_name, chrom_start, chrom_end, sample_list1, sample_list2, size, step, het_threshold)
            if fst_data.size > 0:
                all_fst_data.append(fst_data)
        
        if not all_fst_data:
            return np.array([])
        
        return np.vstack(all_fst_data)

    def get_tajima_d_data(self, chrom_name, chrom_start, chrom_end, size=1, step=1, sample_list=None, het_threshold=None, include_qc=False):
        """
        Get Tajima's D data of target samples using sliding window method

        :param chrom_name: chromosome name
        :param chrom_start: chromosome start position
        :param chrom_end: chromosome end position
        :param size: window size
        :param step: step size
        :param sample_list: a list contained sample names
        :param het_threshold: heterozygosity rate threshold for filtering variants
        :param include_qc: append fixed sample size, retained/excluded SNP counts,
                           and segregating sites (TAJIMA_QC_COLUMNS)
        :return: a numpy array with TAJIMA_COLUMNS, or TAJIMA_QC_COLUMNS

        Keep the selected population fixed and exclude each SNP with any
        missing diploid call in that population. Apply the original Tajima D
        formula to the remaining complete sites; no variable-n approximation.
        """
        windows = self._population_windows(chrom_start, chrom_end, size, step)
        pos, ac, n_chromosomes = self._population_allele_counts(
            chrom_name, chrom_start, chrom_end, sample_list, het_threshold)

        complete = ac.sum(axis=1) == n_chromosomes
        rows = []
        for left, right in windows:
            first, last = np.searchsorted(pos, left), np.searchsorted(pos, right, side='right')
            window_ac = ac[first:last][complete[first:last]]
            counts = int(last - first)
            used = len(window_ac)
            excluded = counts - used
            segregating = int(window_ac.count_segregating())
            value = np.nan
            if n_chromosomes >= 4 and segregating >= 3:
                value = float(allel.tajima_d(window_ac, min_sites=3))
                if not np.isfinite(value):
                    value = np.nan
            rows.append([chrom_name, int(left), int(right), counts, value,
                         n_chromosomes // 2, n_chromosomes, used, excluded,
                         segregating])
        if np.any(~complete):
            logger.warning(f"Tajima's D excludes {np.count_nonzero(~complete)} SNPs with missing "
                           "calls in the selected population; uses complete sites with a fixed "
                           "sample set. Non-random missingness can bias the retained spectrum.")
        if include_qc:
            return np.asarray(rows, dtype=object)
        # Normal API/CLI output uses five columns; extra diagnostics are opt-in.
        return np.asarray([row[:len(TAJIMA_COLUMNS)] for row in rows], dtype=str)

    def get_multi_loci_tajima_d_data(self, loci, size=1, step=1, sample_list=None, het_threshold=None, include_qc=False):
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
            tajima_d_data = self.get_tajima_d_data(chrom_name, chrom_start, chrom_end, size, step, sample_list, het_threshold, include_qc)
            if tajima_d_data.size > 0:
                all_tajima_d_data.append(tajima_d_data)
        
        if not all_tajima_d_data:
            return np.empty((0, len(TAJIMA_QC_COLUMNS) if include_qc else len(TAJIMA_COLUMNS)),
                            dtype=object if include_qc else str)
        
        return np.vstack(all_tajima_d_data)

    def get_ld_data(self, chrom_name, chrom_start, chrom_end):
        """Return a square r² matrix in VCF record order (empty regions: 0 x 0)."""
        logger.info(f"Calculating LD data in target region {chrom_name}:{chrom_start}-{chrom_end}")
        geno_data = self.get_geno_data(chrom_name, chrom_start, chrom_end, 'code1')
        if not geno_data.index.is_unique:
            raise ValueError('Duplicate variant records in LD input')
        return ld_r_squared(geno_data.to_numpy())

    def get_multi_loci_ld_data(self, loci):
        """Return r² for sorted, merged intervals; overlapping requests count once.

        Matrix order follows get_multi_loci_geno_data(normalize_loci(loci)).
        Distinct chromosomes are allowed, but r² is not evidence of physical linkage.
        """
        loci = normalize_loci(loci)
        logger.info(f"Calculating LD data from {len(loci)} regions")
        geno_data = self.get_multi_loci_geno_data(loci, 'code1')
        if not geno_data.index.is_unique:
            raise ValueError('Duplicate variant records in LD input')
        return ld_r_squared(geno_data.to_numpy())

    def get_heterozygosity_data(self, chrom_name, chrom_start, chrom_end):
        """
        Get heterozygosity data of all samples

        :return: a numpy array contained heterozygosity data
        """
        geno_data = self.get_geno_data(chrom_name, chrom_start, chrom_end, 'code2')
        g = allel.GenotypeArray(geno_data.to_numpy().tolist())
        # the rate of observed heterozygosity
        het_o = allel.heterozygosity_observed(g)
        # the expected rate of heterozygosity
        af = g.count_alleles().to_frequencies()
        het_e = allel.heterozygosity_expected(af, ploidy=2)
        return np.column_stack([het_o, het_e])

    def get_inbreeding_coefficient(self, chrom_name, chrom_start, chrom_end):
        """
        Get inbreeding coefficient of all samples

        :return: a numpy array contained inbreeding coefficient
        """
        geno_data = self.get_geno_data(chrom_name, chrom_start, chrom_end, 'code2')
        g = allel.GenotypeArray(geno_data.to_numpy().tolist())
        return allel.inbreeding_coefficient(g)

    def filter_by_heterozygosity(self, geno_data, het_threshold):
        """
        Filter variants by heterozygosity rate
        
        :param geno_data: DataFrame with genotype data (code1 format: 0=homozygous ref, 1=heterozygous, 2=homozygous alt)
        :param het_threshold: heterozygosity rate threshold, variants with het_rate > threshold will be filtered out
        :return: filtered DataFrame
        """
        if het_threshold is None or geno_data.empty:
            return geno_data
        
        logger.info(f"Filtering variants by heterozygosity rate threshold: {het_threshold}")
        
        # Calculate heterozygosity rate for each variant
        het_rates = []
        variants_to_keep = []
        
        for variant_idx in geno_data.index:
            variant_genotypes = geno_data.loc[variant_idx]
            
            # Count heterozygous genotypes (genotype = 1 in code1 format)
            total_samples = len(variant_genotypes)
            missing_samples = (variant_genotypes == -1).sum()  # Missing genotypes
            valid_samples = total_samples - missing_samples
            
            if valid_samples == 0:
                # Skip variants with no valid genotypes
                continue
            
            het_count = (variant_genotypes == 1).sum()
            het_rate = het_count / valid_samples
            
            het_rates.append(het_rate)
            
            # Keep variant if heterozygosity rate <= threshold
            if het_rate <= het_threshold:
                variants_to_keep.append(variant_idx)
        
        # Filter the genotype data
        if len(variants_to_keep) > 0:
            filtered_data = geno_data.loc[variants_to_keep]
        else:
            # Return empty DataFrame with same structure if no variants pass filter
            filtered_data = geno_data.iloc[0:0].copy()
        
        original_count = len(geno_data)
        filtered_count = len(filtered_data)
        removed_count = original_count - filtered_count
        
        logger.info(f"Heterozygosity filtering: {original_count} original variants, "
                   f"{filtered_count} retained, {removed_count} removed")
        
        if len(het_rates) > 0:
            avg_het_rate = np.mean(het_rates)
            max_het_rate = np.max(het_rates)
            min_het_rate = np.min(het_rates)
            logger.info(f"Heterozygosity rates - Average: {avg_het_rate:.4f}, "
                       f"Min: {min_het_rate:.4f}, Max: {max_het_rate:.4f}")
        
        return filtered_data


class MultiPopulationAnalyzer(object):
    """Compare genotype profiles on a common, complete set of VCF records."""

    def __init__(self, gv_list, pop_names):
        self.gv_list, self.pop_names = list(gv_list), list(pop_names)
        if len(self.gv_list) < 2 or len(self.gv_list) != len(self.pop_names):
            raise ValueError("Compare requires at least two VCFs and one name per population")
        if any(not isinstance(name, str) or not name.strip() for name in self.pop_names):
            raise ValueError("Population names must be non-empty strings")
        if len(set(self.pop_names)) != len(self.pop_names):
            raise ValueError("Population names must be unique")
        for name in self.pop_names:
            if (name != name.strip() or name.endswith('.') or name.lower() == 'merged' or
                    any(ch in r'<>:"/\|?*' or ord(ch) < 32 for ch in name)):
                raise ValueError("Population names must be safe filename labels; 'merged' is reserved")

    @staticmethod
    def _merged_loci(loci):
        return normalize_loci(loci)

    @staticmethod
    def _canonical_allele(allele):
        # VCF sequence alleles are case-insensitive; symbolic IDs are not.
        return allele.upper() if allele and set(allele.upper()) <= set('ACGTN') else allele

    def compare_populations(self, loci):
        loci = self._merged_loci(loci)
        pop_geno_data, all_variants, snp_refs = {}, set(), {}
        for i, (gv, name) in enumerate(zip(self.gv_list, self.pop_names), 1):
            data = (gv.get_geno_data(*loci[0]) if len(loci) == 1
                    else gv.get_multi_loci_geno_data(loci)).copy()
            data.index = pd.MultiIndex.from_tuples(
                [(chrom, pos, self._canonical_allele(ref), self._canonical_allele(alt))
                 for chrom, pos, ref, alt in data.index], names=['chrom', 'pos', 'ref', 'alt'])
            if data.index.has_duplicates:
                raise ValueError(f"Duplicate variant records in population {name}")
            if not data.columns.is_unique:
                raise ValueError(f"Duplicate sample names in population {name}")
            for chrom, pos, ref, alt in data.index:
                if ref in ('A', 'C', 'G', 'T') and alt in ('A', 'C', 'G', 'T'):
                    refs = snp_refs.setdefault((chrom, pos), set())
                    refs.add(ref)
                    if len(refs) > 1:
                        raise ValueError(f"Incompatible SNP REF alleles at {chrom}:{pos}; "
                                         "align reference assembly and allele orientation before compare")
            data = data.sort_index()
            data.columns = pd.Index([f"{sample}_{i}" for sample in data.columns], name='samples')
            # Keep empty inputs in presence columns and sample-count summaries.
            pop_geno_data[name] = data
            all_variants.update(data.index)
            logger.info(f"Population {name}: {data.shape[1]} samples, {len(data)} records")

        presence = self._analyze_variant_presence(pop_geno_data, all_variants)
        core = list(presence.index[presence['present_in_all'].to_numpy(dtype=bool)])
        if core:
            combined = pd.concat([data.loc[core] for data in pop_geno_data.values()], axis=1)
            complete = combined.isin([0, 1, 2]).all(axis=1)
            used = list(combined.index[complete])
            excluded = len(core) - len(used)
            if excluded:
                logger.warning(f"Compare excludes {excluded} of {len(core)} shared records with "
                               "missing calls; all populations use the same complete loci")
        else:
            used = []
        logger.info(f"Compare uses {len(used)} complete shared loci for genotype-profile grouping")
        if not used:
            logger.warning("No complete shared loci; no haplotype assignments can be made")

        results = {}
        for name, data in pop_geno_data.items():
            table, groups = self._haplotype_outputs(data.loc[used], f"{name}Hap")
            results[f'{name}_hap_table'], results[f'{name}_hap_groups'] = table, groups
        combined = pd.concat([data.loc[used] for data in pop_geno_data.values()], axis=1)
        results['merged_hap_table'], results['merged_hap_groups'] = self._haplotype_outputs(combined, 'MergedHap')
        results['comparison_summary'] = self._generate_comparison_summary(pop_geno_data, presence)
        results['variant_class'] = presence
        return results

    def _analyze_variant_presence(self, pop_geno_data, all_variants):
        columns = ['variant', 'chrom', 'pos', 'ref', 'alt']
        columns += [f'present_in_{name}' for name in self.pop_names]
        columns += ['present_in_populations', 'num_populations', 'present_in_all',
                    'population_specific', 'specific_to']
        rows = []
        for variant in sorted(all_variants):
            present = [name for name in self.pop_names if variant in pop_geno_data[name].index]
            row = dict(zip(['chrom', 'pos', 'ref', 'alt'], variant))
            row['variant'] = variant
            row.update({f'present_in_{name}': name in present for name in self.pop_names})
            row.update(present_in_populations=','.join(present), num_populations=len(present),
                       present_in_all=len(present) == len(self.pop_names),
                       population_specific=len(present) == 1,
                       specific_to=present[0] if len(present) == 1 else 'None')
            rows.append(row)
        frame = pd.DataFrame(rows, columns=columns).set_index('variant')
        for col in [f'present_in_{name}' for name in self.pop_names] + ['present_in_all', 'population_specific']:
            frame[col] = frame[col].astype(bool)
        frame['num_populations'] = frame['num_populations'].astype(int)
        return frame

    @staticmethod
    def _haplotype_outputs(data, prefix):
        return haplotype_outputs(data, prefix)

    def _generate_comparison_summary(self, pop_geno_data, variant_presence):
        total = len(variant_presence)
        core = int(variant_presence['present_in_all'].sum())
        rows = [
            dict(metric='Total variants', value=total, description='Distinct biallelic records across input VCFs'),
            dict(metric='Core variants', value=core,
                 description='Records in every input VCF, before excluding incomplete loci for grouping'),
            dict(metric='Core variant percentage', value=round(core / total * 100, 2) if total else 0,
                 description='Percentage of distinct input records present in every VCF'),
        ]
        for name, data in pop_geno_data.items():
            rows.extend([
                dict(metric=f'{name} samples', value=len(data.columns), description=f'Input samples in {name}'),
                dict(metric=f'{name} total variants', value=len(data), description=f'Input biallelic records in {name}'),
                dict(metric=f'{name} core variant coverage',
                     value=round(core / len(data) * 100, 2) if len(data) else 0,
                     description=f'Shared input-record percentage in {name}'),
            ])
        return pd.DataFrame(rows, columns=['metric', 'value', 'description'])


class GeneHapNetwork:
    """
    A class for gene haplotype network analysis
    """
    def __init__(self):
        self.edge_data = None
        self.network = None
        self.nodes = []

    @staticmethod
    def _node_key(node):
        return type(node).__name__, str(node)

    def calc_hamming_distance(self, df: pd.DataFrame):
        """
        Calculate the Hamming distance for all pairs of haplotypes

        :param df: DataFrame with haplotypes. each row is a locus. each column is a haplotype. values are 0, 1 or 2.
        :return: a list of tuples, each tuple contains the source, target and weight of an edge
        """
        if not isinstance(df, pd.DataFrame) or df.shape[1] == 0:
            raise ValueError("Haplotype table must contain at least one haplotype")
        if not df.columns.is_unique or df.columns.hasnans:
            raise ValueError("Haplotype names must be unique and non-missing")
        if any(not str(name).strip() or any(c in str(name) for c in '\t\r\n') for name in df.columns):
            raise ValueError("Haplotype names must be non-empty and contain no tab/newline")
        haplotypes = sorted(df.columns, key=self._node_key)
        try:
            values = df.loc[:, haplotypes].apply(pd.to_numeric, errors='raise').to_numpy(dtype=float, na_value=np.nan)
        except (TypeError, ValueError) as exc:
            raise ValueError("Haplotype states must be 0, 1, 2, -1 or missing") from exc
        if not np.all(np.isnan(values) | np.isin(values, [-1, 0, 1, 2])):
            raise ValueError("Haplotype states must be 0, 1, 2, -1 or missing")
        # All pairs must use the same loci, preserving an ordinary Hamming
        # metric. Neither -1 nor NaN is a biological genotype state.
        complete = np.all(np.isfinite(values) & (values != -1), axis=1)
        if not np.any(complete):
            raise ValueError("No loci with complete calls across all haplotypes")
        excluded = int(np.count_nonzero(~complete))
        if excluded:
            logger.warning(f"Network distances exclude {excluded} of {len(values)} loci "
                           "with missing calls across haplotypes")
        values = values[complete].T
        return [(haplotypes[i], haplotypes[j], int(np.count_nonzero(values[i] != values[j])))
                for i in range(len(haplotypes)) for j in range(i + 1, len(haplotypes))]

    def add_edge(self, df: pd.DataFrame):
        """
        Add an edge to the network

        :param df: DataFrame with haplotypes. each row is a locus. each column is a haplotype. values are 0, 1 or 2.
        """
        edges = self.calc_hamming_distance(df)
        self.nodes = sorted(df.columns, key=self._node_key)
        self.edge_data = edges
        self.network = None

    def _distance_graph(self):
        """Build a graph in canonical order, including nodes with no edges."""
        from numbers import Real
        if self.edge_data is None:
            raise ValueError("No edge data available. Call add_edge() first.")
        edges, nodes = {}, set(self.nodes)
        for u, v, weight in self.edge_data:
            if u == v:
                raise ValueError("Distance edges must connect distinct nodes")
            if not isinstance(weight, Real) or not np.isfinite(weight) or weight < 0:
                raise ValueError("Distance weights must be finite non-negative numbers")
            u, v = sorted((u, v), key=self._node_key)
            if (u, v) in edges and edges[u, v] != weight:
                raise ValueError("Conflicting duplicate distance edges")
            edges[u, v] = weight
            nodes.update((u, v))
        graph = nx.Graph()
        graph.add_nodes_from(sorted(nodes, key=self._node_key))
        ordered = sorted(edges.items(), key=lambda item: (item[1], self._node_key(item[0][0]), self._node_key(item[0][1])))
        graph.add_weighted_edges_from((u, v, weight) for (u, v), weight in ordered)
        return graph

    def minimum_spanning_tree(self):
        """
        Construct the minimum spanning tree from distance matrix
        """
        return nx.minimum_spanning_tree(self._distance_graph(), algorithm='kruskal', weight='weight')

    def minimum_spanning_network(self):
        """
        Return the union of all minimum spanning forests, retaining tied edges.

        Process equal-weight edges as a batch. An edge is eligible exactly
        when its endpoints are disconnected using strictly lighter edges.
        Union components only after every edge in the batch has been tested.
        """
        from itertools import groupby
        graph = self._distance_graph()
        network = nx.Graph()
        network.add_nodes_from(graph.nodes())
        components = nx.utils.UnionFind(graph.nodes())
        edges = sorted(graph.edges(data=True), key=lambda edge: edge[2]['weight'])
        for _, group in groupby(edges, key=lambda edge: edge[2]['weight']):
            batch = list(group)
            for u, v, data in batch:
                if components[u] != components[v]:
                    network.add_edge(u, v, weight=data['weight'])
            for u, v, _ in batch:
                components.union(u, v)
        logger.info(f"MSN contains {network.number_of_nodes()} nodes and "
                    f"{network.number_of_edges()} edges, including equal-weight alternatives")
        return network

    def create_network(self, method: str = "MST"):
        """
        Create the network from distance matrix using different methods

        :param method: Method to construct the network. Available methods:
                      - "MST": Minimum Spanning Tree
                      - "MSN": Minimum Spanning Network
        """
        if self.edge_data is None:
            raise ValueError("No edge data available. Call add_edge() first.")
        
        logger.info(f"Creating network using method: {method}")
        
        if method.upper() == "MST":
            self.network = self.minimum_spanning_tree()
        elif method.upper() == "MSN":
            self.network = self.minimum_spanning_network()
        else:
            available_methods = ["MST", "MSN"]
            raise ValueError(f"Invalid method: {method}. Available methods: {available_methods}")

    def save_network(self, file_path: str, format: str = 'auto'):
        """
        Save the network to a file in various formats

        :param file_path: The path to save the network
        :param format: File format ('gexf', 'gml', 'edgelist', 'csv', 'auto'). 
                      If 'auto', format is determined by file extension.
        """
        from pathlib import Path
        
        if self.network is None:
            raise ValueError("No network created yet. Please call create_network() first.")
        
        # Determine format from file extension if auto
        if format == 'auto':
            ext = Path(file_path).suffix.lower()
            format_map = {
                '.gexf': 'gexf',
                '.gml': 'gml', 
                '.txt': 'edgelist',
                '.edgelist': 'edgelist',
                '.csv': 'csv'
            }
            format = format_map.get(ext, 'csv')
        
        logger.info(f"Saving network in {format.upper()} format to: {file_path}")
        
        # Save in requested format
        if format == 'gexf':
            # Add node attributes for better visualization
            for node in self.network.nodes():
                self.network.nodes[node]['label'] = node
                self.network.nodes[node]['id'] = node
            
            # Add edge attributes
            for u, v, data in self.network.edges(data=True):
                data['label'] = f"weight: {data['weight']}"
                data['type'] = 'undirected'
            
            nx.write_gexf(self.network, file_path)
            
        elif format == 'gml':
            # GML format - ensure all attributes are valid
            network_copy = self.network.copy()
            
            # Add node attributes
            for node in network_copy.nodes():
                network_copy.nodes[node]['label'] = str(node)
                network_copy.nodes[node]['id'] = str(node)
            
            # Ensure edge weights are numeric
            for u, v, data in network_copy.edges(data=True):
                if 'weight' in data:
                    data['weight'] = float(data['weight'])
            
            nx.write_gml(network_copy, file_path)
            
        elif format == 'edgelist':
            # Simple edge list format
            import json
            with open(file_path, 'w', encoding='utf-8') as f:
                for node in nx.isolates(self.network):
                    # A comment keeps normal edge-list readers compatible;
                    # hastat's reader restores the otherwise unrepresentable node.
                    f.write('#node\t' + json.dumps(str(node), ensure_ascii=False) + '\n')
                for u, v, data in self.network.edges(data=True):
                    weight = data.get('weight', 1)
                    f.write(f"{u}\t{v}\t{weight}\n")
                    
        elif format == 'csv':
            # Original CSV format with additional information
            results = []
            
            # Save edge information from the constructed network
            for u, v, data in self.network.edges(data=True):
                results.append({
                    'source': u,
                    'target': v, 
                    'weight': data.get('weight', 0),
                    'type': 'edge'
                })

            for node in nx.isolates(self.network):
                results.append({'source': node, 'target': '', 'weight': '', 'type': 'node'})
            # Existing columns also represent isolated nodes, without self-loops.
            results_df = pd.DataFrame(results, columns=['source', 'target', 'weight', 'type'])
            results_df.to_csv(file_path, index=False)
            
        else:
            raise ValueError(f"Unsupported format: {format}. Supported formats: gexf, gml, edgelist, csv")
        
        logger.info(f"Network saved successfully with {self.network.number_of_nodes()} nodes and {self.network.number_of_edges()} edges")
