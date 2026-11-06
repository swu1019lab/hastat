"""Variant extraction, filtering and genomic intervals."""
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


class VariantData:
    """Read diploid genotype data from an indexed VCF handle."""

    def __init__(self, vcf_in) -> None:
        """
        :param vcf_in: VariantFile object
        """
        self.vcf_in = vcf_in

    def extract_genotypes(self, chrom_name, chrom_start, chrom_end, code='code1', het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
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

    def extract_multi_loci_genotypes(self, loci, code='code1', het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
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
            region_data = self.extract_genotypes(chrom_name, chrom_start, chrom_end, code, het_threshold, maf_threshold, max_missing_threshold)
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
