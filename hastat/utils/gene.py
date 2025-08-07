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
from scipy.spatial.distance import squareform
from hastat.hastat import logger


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
        if gene.strand == "+":
            start = gene.start - upstream if upstream < gene.start else 0
            end = gene.end + downstream
        else:
            start = gene.start - downstream if downstream < gene.start else 0
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

    def hap_table(self, chrom_name, chrom_start, chrom_end, het_threshold=None):
        """
        Get haplotypes table of a gene region

        :param chrom_name: chromosome name
        :param chrom_start: chromosome start position
        :param chrom_end: chromosome end position
        :param het_threshold: heterozygosity rate threshold for filtering variants
        :return: a DataFrame contained haplotypes table
        """
        logger.info("Initializing the GeneVariant object")

        df = self.get_geno_data(chrom_name, chrom_start, chrom_end, het_threshold=het_threshold)
        # get unique haplotypes and their counts
        hap, counts = np.unique(df.T.values, axis=0, return_counts=True)

        # Sort haplotypes by sample count (descending order)
        # Hap1 should have the most samples, HapN should have the least
        sorted_indices = np.argsort(counts)[::-1]  # Sort indices in descending order of counts
        hap_sorted = hap[sorted_indices]

        # haplotypes table, including name, genotype code and size
        hap_table = pd.DataFrame(hap_sorted, columns=df.index)
        # hap_table['size'] = _  # Could be used to add sample counts to table
        hap_table.index = "Hap" + (hap_table.index + 1).astype(str)
        # hap_table.index.name = 'haplotypes'
        hap_table = hap_table.T.reset_index()

        # code number should transform into code character
        # i: chr, pos, ref, alt
        # for i, s in hap_table.iloc[:, :-1].items():
        #     ref, alt = i[2], i[3]
        #     hap_table.loc[:, i] = s.map({0: ref + ref, 1: ref + alt, 2: alt + alt})
        # hap_table.columns = hap_table.columns.droplevel(['ref', 'alt'])
        return hap_table

    def multi_loci_hap_table(self, loci, het_threshold=None):
        """
        Get haplotypes table from multiple regions (for homologous genes)

        :param loci: list of tuples containing (chrom_name, chrom_start, chrom_end)
        :param het_threshold: heterozygosity rate threshold for filtering variants
        :return: a DataFrame contained haplotypes table
        """
        logger.info(f"Processing haplotype table for {len(loci)} regions")
        
        df = self.get_multi_loci_geno_data(loci, het_threshold=het_threshold)
        # get unique haplotypes and their counts
        hap, counts = np.unique(df.T.values, axis=0, return_counts=True)

        # Sort haplotypes by sample count (descending order)
        # Hap1 should have the most samples, HapN should have the least
        sorted_indices = np.argsort(counts)[::-1]  # Sort indices in descending order of counts
        hap_sorted = hap[sorted_indices]

        # haplotypes table, including name, genotype code and size
        hap_table = pd.DataFrame(hap_sorted, columns=df.index)
        hap_table.index = "Hap" + (hap_table.index + 1).astype(str)
        hap_table = hap_table.T.reset_index()
        
        return hap_table

    def hap_groups(self, chrom_name, chrom_start, chrom_end, het_threshold=None):
        """
        Get haplotypes group of a gene region

        :param chrom_name: chromosome name
        :param chrom_start: chromosome start position
        :param chrom_end: chromosome end position
        :param het_threshold: heterozygosity rate threshold for filtering variants
        :return: a DataFrame contained haplotypes group
        """
        logger.info("Initializing the GeneVariant object")

        df = self.get_geno_data(chrom_name, chrom_start, chrom_end, het_threshold=het_threshold)
        # get unique haplotypes, their counts and sample indices
        _, indices, counts = np.unique(df.T.values, axis=0, return_inverse=True, return_counts=True)

        # Sort haplotypes by sample count (descending order)
        # Hap1 should have the most samples, HapN should have the least
        sorted_indices = np.argsort(counts)[::-1]  # Sort indices in descending order of counts
        
        # Create mapping from haplotype indices to their ranks (sorted by count)
        count_to_rank = np.zeros(len(sorted_indices), dtype=int)
        count_to_rank[sorted_indices] = np.arange(len(sorted_indices))
        
        # Remap the sample indices according to the new sorted order
        remapped_indices = count_to_rank[indices]

        # haplotypes group for each sample
        # add haplotype labels, e.g. Hap1, Hap2, ..., HapN
        hap_groups = pd.Series(remapped_indices, index=df.columns)
        hap_groups = "Hap" + (hap_groups + 1).astype(str)
        hap_groups.name = 'haplotypes'
        return hap_groups.reset_index()

    def multi_loci_hap_groups(self, loci, het_threshold=None):
        """
        Get haplotypes group from multiple regions (for homologous genes)

        :param loci: list of tuples containing (chrom_name, chrom_start, chrom_end)
        :param het_threshold: heterozygosity rate threshold for filtering variants
        :return: a DataFrame contained haplotypes group
        """
        logger.info(f"Processing haplotype groups for {len(loci)} regions")
        
        df = self.get_multi_loci_geno_data(loci, het_threshold=het_threshold)
        # get unique haplotypes, their counts and sample indices
        _, indices, counts = np.unique(df.T.values, axis=0, return_inverse=True, return_counts=True)

        # Sort haplotypes by sample count (descending order)
        # Hap1 should have the most samples, HapN should have the least
        sorted_indices = np.argsort(counts)[::-1]  # Sort indices in descending order of counts
        
        # Replace the sample indices with their ranks (sorted by count)
        count_to_rank = np.zeros(len(sorted_indices), dtype=int)
        count_to_rank[sorted_indices] = np.arange(len(sorted_indices))
        
        # Remap the sample indices according to the new sorted order
        remapped_indices = count_to_rank[indices]

        # haplotypes group for each sample
        hap_groups = pd.Series(remapped_indices, index=df.columns)
        hap_groups = "Hap" + (hap_groups + 1).astype(str)
        hap_groups.name = 'haplotypes'
        return hap_groups.reset_index()

    def calc_hap_freq(self, chrom_name, chrom_start, chrom_end, sample_groups=None, het_threshold=None):
        """
        Calculate haplotypes frequency of a gene region

        :param chrom_name: chromosome name
        :param chrom_start: chromosome start position
        :param chrom_end: chromosome end position
        :param sample_groups: A sample group Dataframe with two columns only: samples and groups
        :param het_threshold: heterozygosity rate threshold for filtering variants
        :return: a DataFrame contained haplotypes frequency
        """
        hap_groups = self.hap_groups(chrom_name, chrom_start, chrom_end, het_threshold=het_threshold)
        if sample_groups is not None:
            return hap_groups.merge(sample_groups, how='inner', on='samples').groupby(
                ['haplotypes', 'groups']).count().unstack(fill_value=0)
        else:
            return hap_groups.groupby('haplotypes').count()

    def calc_multi_loci_hap_freq(self, loci, sample_groups=None, het_threshold=None):
        """
        Calculate haplotypes frequency from multiple regions (for homologous genes)

        :param loci: list of tuples containing (chrom_name, chrom_start, chrom_end)
        :param sample_groups: A sample group Dataframe with two columns only: samples and groups
        :param het_threshold: heterozygosity rate threshold for filtering variants
        :return: a DataFrame contained haplotypes frequency
        """
        hap_groups = self.multi_loci_hap_groups(loci, het_threshold=het_threshold)
        if sample_groups is not None:
            return hap_groups.merge(sample_groups, how='inner', on='samples').groupby(
                ['haplotypes', 'groups']).count().unstack(fill_value=0)
        else:
            return hap_groups.groupby('haplotypes').count()

    def get_geno_data(self, chrom_name, chrom_start, chrom_end, code='code1', het_threshold=None):
        """
        Get genotype data within a gene region

        :param chrom_name: chromosome name
        :param chrom_start: chromosome start position
        :param chrom_end: chromosome end position
        :param code: genotype code, code1 or code2, default code1, code1 is presented as 0/1/2,
                     code2 is presented as [0, 0], [0, 1] and [1, 1]
        :param het_threshold: heterozygosity rate threshold for filtering variants
        :return: a DataFrame contained genotype data
        """
        allele_data_set = {'code1': defaultdict(list), 'code2': defaultdict(list)}
        allele_data_set['code1']['index_names'] = ['chrom', 'pos', 'ref', 'alt']
        allele_data_set['code1']['column_names'] = ['samples']
        allele_data_set['code1']['columns'] = list(self.vcf_in.header.samples)
        allele_data_set['code2']['index_names'] = ['chrom', 'pos', 'ref', 'alt']
        allele_data_set['code2']['column_names'] = ['samples']
        allele_data_set['code2']['columns'] = list(self.vcf_in.header.samples)
        # Access snp/InDel loci
        # 0-based, half-open
        try:
            _ = next(self.vcf_in.fetch(chrom_name, chrom_start, chrom_end))
        except StopIteration:
            logger.warning("No SNP or InDel loci existed in target region!")
        except ValueError as e:
            logger.error("{}".format(e))
            logger.error("\tAn error occurred in {}:{}-{}".format(chrom_name, chrom_start, chrom_end))
            raise e
        else:
            for rec in self.vcf_in.fetch(chrom_name, chrom_start, chrom_end):
                # print(rec.chrom, rec.pos, rec.alleles)
                # only keep bi-allelic SNP/InDel loci
                if len(rec.alleles) != 2:
                    continue
                # genotype code was represented by number, e.g. 0, 1, and 2
                # Missing genotypes are represented by -1
                geno_code1 = [sum(s.allele_indices) if isinstance(s.allele_indices[0], int) else -1 for s in
                              rec.samples.values()]
                allele_data_set['code1']['index'].append((rec.chrom, rec.pos, rec.ref, rec.alts[0]))
                allele_data_set['code1']['data'].append(geno_code1)
                # genotype code was represented by paired number, e.g., [0, 0], [0, 1] and [1, 1]
                # Missing genotypes are represented by [-1, -1]
                geno_code2 = [list(s.allele_indices) if isinstance(s.allele_indices[0], int) else [-1, -1] for s in
                              rec.samples.values()]
                allele_data_set['code2']['index'].append((rec.chrom, rec.pos, rec.ref, rec.alts[0]))
                allele_data_set['code2']['data'].append(geno_code2)
        
        geno_data = pd.DataFrame.from_dict(allele_data_set[code], orient='tight')
        
        # Apply heterozygosity filtering if threshold is provided
        # For code2 format, we need to convert to code1, filter, then extract the same variants
        if het_threshold is not None and not geno_data.empty:
            if code == 'code1':
                geno_data = self.filter_by_heterozygosity(geno_data, het_threshold)
            elif code == 'code2':
                # Get code1 data for filtering
                code1_data = pd.DataFrame.from_dict(allele_data_set['code1'], orient='tight')
                if not code1_data.empty:
                    filtered_code1_data = self.filter_by_heterozygosity(code1_data, het_threshold)
                    # Keep only the variants that passed filtering in code1
                    if not filtered_code1_data.empty:
                        geno_data = geno_data.loc[filtered_code1_data.index]
                    else:
                        # Return empty DataFrame with same structure if no variants pass filter
                        geno_data = geno_data.iloc[0:0].copy()
        
        return geno_data

    def get_multi_loci_geno_data(self, loci, code='code1', het_threshold=None):
        """
        Get genotype data from multiple regions (for homologous genes)

        :param loci: list of tuples containing (chrom_name, chrom_start, chrom_end)
        :param code: genotype code, code1 or code2, default code1
        :param het_threshold: heterozygosity rate threshold for filtering variants
        :return: a DataFrame contained genotype data from all regions
        """
        logger.info(f"Extracting genotype data from {len(loci)} regions")
        
        all_data = []
        for chrom_name, chrom_start, chrom_end in loci:
            logger.info(f"Processing region: {chrom_name}:{chrom_start}-{chrom_end}")
            region_data = self.get_geno_data(chrom_name, chrom_start, chrom_end, code)
            if not region_data.empty:
                # Apply heterozygosity filtering if threshold is provided
                if het_threshold is not None and code == 'code1':
                    region_data = self.filter_by_heterozygosity(region_data, het_threshold)
                all_data.append(region_data)
        
        if not all_data:
            logger.warning("No genotype data found in any of the target regions")
            return pd.DataFrame()
        
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
        logger.info(f"Calculating Pi data of {len(sample_list)} samples using sliding window method")
        geno_data = self.get_geno_data(chrom_name, chrom_start, chrom_end, het_threshold=het_threshold)
        if sample_list is not None:
            data = geno_data.loc[:, geno_data.columns.isin(sample_list)].to_numpy()
        else:
            data = geno_data.to_numpy()
        h = allel.HaplotypeArray(data)
        ac = h.count_alleles()
        pos = geno_data.index.get_level_values('pos').to_list()
        pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac, size, step=step)
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
        logger.info(
            f"Calculating Fst data between {len(sample_list1)} and {len(sample_list2)} samples "
            f"using sliding window method")
        geno_data = self.get_geno_data(chrom_name, chrom_start, chrom_end, 'code2', het_threshold=het_threshold)
        # Returns -1 for unmatched values
        idx = geno_data.columns
        index_group1 = idx.get_indexer(idx.where(idx.isin(sample_list1)).dropna())
        index_group2 = idx.get_indexer(idx.where(idx.isin(sample_list2)).dropna())
        sub_pops = [index_group1, index_group2]
        pos = geno_data.index.get_level_values('pos').to_list()
        g = allel.GenotypeArray(geno_data.to_numpy().tolist())
        _ = np.seterr(divide='ignore')
        # fst: shape (n_windows,), windows: shape (n_windows, 2), counts: shape (n_windows,)
        fst, windows, counts = allel.windowed_weir_cockerham_fst(pos, g, sub_pops, size, step=step)
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

    def get_ld_data(self, chrom_name, chrom_start, chrom_end):
        """
        Get LD data of all samples

        :return: a numpy array contained LD data
        """
        logger.info("Calculating LD data in target region {chrom_name}:{chrom_start}-{chrom_end}")
        geno_data = self.get_geno_data(chrom_name, chrom_start, chrom_end, 'code1')
        r = allel.rogers_huff_r(geno_data.to_numpy())
        return squareform(np.power(r, 2))

    def get_multi_loci_ld_data(self, loci):
        """
        Get LD data from multiple regions (for homologous genes)

        :param loci: list of tuples containing (chrom_name, chrom_start, chrom_end)
        :return: a numpy array contained LD data
        """
        logger.info(f"Calculating LD data from {len(loci)} regions")
        
        geno_data = self.get_multi_loci_geno_data(loci, 'code1')
        if geno_data.empty:
            return np.array([])
        
        r = allel.rogers_huff_r(geno_data.to_numpy())
        return squareform(np.power(r, 2))

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
    """
    A class for multi-population haplotype comparison analysis
    """
    
    def __init__(self, gv_list, pop_names):
        """
        Initialize MultiPopulationAnalyzer
        
        :param gv_list: list of GeneVariant objects
        :param pop_names: list of population names
        """
        self.gv_list = gv_list
        self.pop_names = pop_names
        logger.info(f"Initialized MultiPopulationAnalyzer with {len(gv_list)} populations")
    
    def compare_populations(self, loci):
        """
        Perform multi-population haplotype comparison analysis
        
        :param loci: list of tuples containing (chrom_name, chrom_start, chrom_end)
        :return: dictionary containing analysis results
        """
        logger.info("Starting multi-population haplotype comparison analysis")
        
        # Extract genotype data from each population
        pop_geno_data = {}
        all_variants = set()
        
        for i, (gv, pop_name) in enumerate(zip(self.gv_list, self.pop_names)):
            logger.info(f"Processing population {pop_name}")
            
            # Get genotype data for current population
            if len(loci) == 1:
                geno_data = gv.get_geno_data(*loci[0])
            else:
                geno_data = gv.get_multi_loci_geno_data(loci)
            
            if geno_data.empty:
                logger.warning(f"No genotype data found for population {pop_name}")
                continue
            
            # Add population suffix to sample names
            geno_data.columns = [f"{col}_{i+1}" for col in geno_data.columns]
            
            pop_geno_data[pop_name] = geno_data
            all_variants.update(geno_data.index)
            
            logger.info(f"Population {pop_name}: {geno_data.shape[1]} samples, {geno_data.shape[0]} variants")
        
        # Identify core and specific variants
        variant_presence = self._analyze_variant_presence(pop_geno_data, all_variants)
        
        # Generate results
        results = {}
        
        # Core sites analysis (variants present in all populations)
        core_variants = variant_presence[variant_presence['present_in_all']].index
        if len(core_variants) > 0:
            logger.info(f"Found {len(core_variants)} core variants shared across all populations")
            
            # Haplotype analysis based on core variants for each population
            core_results = self._analyze_population_core_haplotype(pop_geno_data, core_variants)
            results.update(core_results)
            
            # Haplotype analysis based on core variants for merged populations
            merged_results = self._analyze_merged_populations_core_haplotype(pop_geno_data, core_variants)
            results.update(merged_results)
        else:
            logger.warning("No core variants found across all populations")
            # Add empty results for each population
            for pop_name in self.pop_names:
                results[f'{pop_name}_hap_table'] = pd.DataFrame()
                results[f'{pop_name}_hap_groups'] = pd.DataFrame()
            # Add empty results for merged populations
            results['merged_hap_table'] = pd.DataFrame()
            results['merged_hap_groups'] = pd.DataFrame()

        # Summary and classification results
        results['comparison_summary'] = self._generate_comparison_summary(pop_geno_data, variant_presence)
        results['variant_class'] = variant_presence
        
        logger.info("Multi-population comparison analysis completed")
        return results
    
    def _analyze_variant_presence(self, pop_geno_data, all_variants):
        """
        Analyze which variants are present in which populations
        """
        variant_presence = []
        
        for variant in all_variants:
            presence_info = {
                'variant': variant,
                'chrom': variant[0],
                'pos': variant[1],
                'ref': variant[2],
                'alt': variant[3]
            }
            
            present_in_pops = []
            for pop_name, geno_data in pop_geno_data.items():
                if variant in geno_data.index:
                    present_in_pops.append(pop_name)
                    presence_info[f'present_in_{pop_name}'] = True
                else:
                    presence_info[f'present_in_{pop_name}'] = False
            
            presence_info['present_in_populations'] = ','.join(present_in_pops)
            presence_info['num_populations'] = len(present_in_pops)
            presence_info['present_in_all'] = len(present_in_pops) == len(self.pop_names)
            presence_info['population_specific'] = len(present_in_pops) == 1
            
            if presence_info['population_specific']:
                presence_info['specific_to'] = present_in_pops[0]
            else:
                presence_info['specific_to'] = 'None'
            
            variant_presence.append(presence_info)
        
        return pd.DataFrame(variant_presence).set_index('variant')
    
    def _analyze_merged_populations_core_haplotype(self, pop_geno_data, core_variants):
        """
        Analyze haplotypes for all merged populations using core variants only
        """
        logger.info(f"Analyzing merged populations with {len(core_variants)} core variants")
        
        # Combine core variant data from all populations
        all_samples = []
        core_data_list = []
        
        for pop_name, geno_data in pop_geno_data.items():
            core_data = geno_data.loc[core_variants]
            core_data_list.append(core_data)
            all_samples.extend(geno_data.columns)
        
        # Combine all core data
        combined_core_data = pd.concat(core_data_list, axis=1)
        
        # Generate haplotype table and groups for merged populations based on core variants
        hap, indices, counts = np.unique(combined_core_data.T.values, axis=0, return_inverse=True, return_counts=True)
        
        # Sort haplotypes by sample count (descending order)
        # MergedHap1 should have the most samples, MergedHapN should have the least
        sorted_indices = np.argsort(counts)[::-1]  # Sort indices in descending order of counts
        hap_sorted = hap[sorted_indices]
        
        # Create mapping from haplotype indices to their ranks (sorted by count)
        count_to_rank = np.zeros(len(sorted_indices), dtype=int)
        count_to_rank[sorted_indices] = np.arange(len(sorted_indices))
        
        # Remap the sample indices according to the new sorted order
        remapped_indices = count_to_rank[indices]
        
        # Haplotype table
        merged_hap_table = pd.DataFrame(hap_sorted, columns=combined_core_data.index)
        merged_hap_table.index = "MergedHap" + (merged_hap_table.index + 1).astype(str)
        merged_hap_table = merged_hap_table.T.reset_index()
        
        # Haplotype groups
        merged_hap_groups = pd.Series(remapped_indices, index=combined_core_data.columns)
        merged_hap_groups = "MergedHap" + (merged_hap_groups + 1).astype(str)
        merged_hap_groups.name = 'haplotypes'
        merged_hap_groups = merged_hap_groups.reset_index()

        logger.info(f"Merged populations: {len(hap)} unique haplotypes from {len(combined_core_data.columns)} samples")
        
        return {
            'merged_hap_table': merged_hap_table,
            'merged_hap_groups': merged_hap_groups
        }
    
    def _analyze_population_core_haplotype(self, pop_geno_data, core_variants):
        """
        Analyze haplotypes within each population separately using core variants only
        """
        results = {}
        
        for pop_name, geno_data in pop_geno_data.items():
            logger.info(f"Analyzing {pop_name} haplotypes based on core variants")
            
            if geno_data.empty:
                continue
            
            # Extract core variants data for this population
            core_data = geno_data.loc[core_variants]
            
            # Generate haplotype table and groups for this population using core variants
            hap, indices, counts = np.unique(core_data.T.values, axis=0, return_inverse=True, return_counts=True)
            
            # Sort haplotypes by sample count (descending order)
            # PopHap1 should have the most samples, PopHapN should have the least
            sorted_indices = np.argsort(counts)[::-1]  # Sort indices in descending order of counts
            hap_sorted = hap[sorted_indices]
            
            # Create mapping from haplotype indices to their ranks (sorted by count)
            count_to_rank = np.zeros(len(sorted_indices), dtype=int)
            count_to_rank[sorted_indices] = np.arange(len(sorted_indices))
            
            # Remap the sample indices according to the new sorted order
            remapped_indices = count_to_rank[indices]
            
            # Haplotype table
            pop_hap_table = pd.DataFrame(hap_sorted, columns=core_data.index)
            pop_hap_table.index = f"{pop_name}Hap" + (pop_hap_table.index + 1).astype(str)
            pop_hap_table = pop_hap_table.T.reset_index()
            
            # Haplotype groups
            pop_hap_groups = pd.Series(remapped_indices, index=core_data.columns)
            pop_hap_groups = f"{pop_name}Hap" + (pop_hap_groups + 1).astype(str)
            pop_hap_groups.name = 'haplotypes'
            pop_hap_groups = pop_hap_groups.reset_index()
            
            results[f'{pop_name}_hap_table'] = pop_hap_table
            results[f'{pop_name}_hap_groups'] = pop_hap_groups
            
            logger.info(f"Population {pop_name}: {len(np.unique(hap, axis=0))} unique haplotypes from {core_data.shape[1]} samples")
        
        return results
    
    def _generate_comparison_summary(self, pop_geno_data, variant_presence):
        """
        Generate summary statistics for the comparison (focusing on core variants)
        """
        summary_data = []
        
        # Overall statistics
        total_variants = len(variant_presence)
        core_variants = len(variant_presence[variant_presence['present_in_all']])
        
        summary_data.append({
            'metric': 'Total variants',
            'value': total_variants,
            'description': 'Total number of variants across all populations'
        })
        
        summary_data.append({
            'metric': 'Core variants',
            'value': core_variants,
            'description': 'Variants present in all populations (used for haplotype analysis)'
        })
        
        summary_data.append({
            'metric': 'Core variant percentage',
            'value': round((core_variants / total_variants * 100), 2) if total_variants > 0 else 0,
            'description': 'Percentage of variants that are core variants'
        })
        
        # Population statistics
        for pop_name, geno_data in pop_geno_data.items():
            if geno_data.empty:
                continue
                
            pop_variants = len(geno_data.index)
            pop_samples = len(geno_data.columns)
            
            summary_data.extend([
                {
                    'metric': f'{pop_name} samples',
                    'value': pop_samples,
                    'description': f'Number of samples in population {pop_name}'
                },
                {
                    'metric': f'{pop_name} total variants',
                    'value': pop_variants,
                    'description': f'Total number of variants in population {pop_name}'
                },
                {
                    'metric': f'{pop_name} core variant coverage',
                    'value': round((core_variants / pop_variants * 100), 2) if pop_variants > 0 else 0,
                    'description': f'Percentage of core variants relative to total variants in {pop_name}'
                }
            ])
        
        return pd.DataFrame(summary_data)


class GeneHapNetwork:
    """
    A class for gene haplotype network analysis
    """
    def __init__(self):
        self.edge_data = None
        self.network = None

    def calc_hamming_distance(self, df: pd.DataFrame):
        """
        Calculate the Hamming distance for all pairs of haplotypes

        :param df: DataFrame with haplotypes. each row is a locus. each column is a haplotype. values are 0, 1 or 2.
        :return: a list of tuples, each tuple contains the source, target and weight of an edge
        """
        # Get all haplotype names
        haplotypes = df.columns.tolist()
        # Transpose dataframe so each row represents a haplotype
        df_t = df.T
        
        # Initialize edge data list
        edges = []
        
        # Calculate Hamming distance between all haplotype pairs
        for i in range(len(haplotypes)):
            for j in range(i+1, len(haplotypes)):
                # Calculate Hamming distance between two haplotypes
                distance = np.sum(df_t.iloc[i] != df_t.iloc[j])
                edges.append((haplotypes[i], haplotypes[j], distance))
        
        return edges

    def add_edge(self, df: pd.DataFrame):
        """
        Add an edge to the network

        :param df: DataFrame with haplotypes. each row is a locus. each column is a haplotype. values are 0, 1 or 2.
        """
        self.edge_data = self.calc_hamming_distance(df)

    def minimum_spanning_tree(self):
        """
        Construct the minimum spanning tree from distance matrix
        """
        G = nx.Graph()
        G.add_weighted_edges_from(self.edge_data)
        T = nx.minimum_spanning_tree(G)
        return T

    def minimum_spanning_network(self, threshold_factor=1.2):
        """
        Construct the minimum spanning network from distance matrix
        
        :param threshold_factor: Threshold factor for adding edges to MSN. 
                               Edges with weight <= path_weight * threshold_factor will be added.
        """
        # Use existing method to construct the minimum spanning tree
        MST = self.minimum_spanning_tree()
        
        # Check for multiple connected components
        num_components = nx.number_connected_components(MST)
        if num_components > 1:
            logger.info(f"Network has {num_components} disconnected components. "
                       f"MSN will preserve this natural separation.")
        
        # Create MSN as a copy of MST
        MSN = MST.copy()
        
        # Get all edge weights in MST for quick lookup
        mst_edges = set()
        for u, v, d in MST.edges(data=True):
            mst_edges.add((min(u, v), max(u, v)))
        
        # Sort all edges by weight for MSN construction
        all_edges = sorted(self.edge_data, key=lambda x: x[2])
        
        # For each edge not in MST, check if it should be added to MSN
        for u, v, weight in all_edges:
            edge = (min(u, v), max(u, v))
            
            # Skip if edge is already in MST
            if edge in mst_edges:
                continue
                
            # Find the shortest path weight between u and v in MST
            try:
                mst_path_weight = nx.shortest_path_length(MST, u, v, weight='weight')
                
                # If the direct edge weight is not much larger than MST path weight, add it
                if weight <= mst_path_weight * threshold_factor:
                    MSN.add_edge(u, v, weight=weight)
            except nx.NetworkXNoPath:
                # If no path exists in MST, nodes are in different connected components
                # Do not force connection - preserve the natural structure
                continue
        
        logger.info(f"MSN construction: MST has {MST.number_of_edges()} edges, "
                   f"MSN has {MSN.number_of_edges()} edges (added {MSN.number_of_edges() - MST.number_of_edges()} edges)")
        
        return MSN

    def create_network(self, method: str = "MST", **kwargs):
        """
        Create the network from distance matrix using different methods

        :param method: Method to construct the network. Available methods:
                      - "MST": Minimum Spanning Tree
                      - "MSN": Minimum Spanning Network
        :param kwargs: Additional parameters for specific methods:
                      - threshold_factor: for MSN (default: 1.2)
        """
        if self.edge_data is None:
            raise ValueError("No edge data available. Call add_edge() first.")
        
        logger.info(f"Creating network using method: {method}")
        
        if method.upper() == "MST":
            self.network = self.minimum_spanning_tree()
        elif method.upper() == "MSN":
            threshold_factor = kwargs.get('threshold_factor', 1.2)
            self.network = self.minimum_spanning_network(threshold_factor)
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
            logger.error("No network created yet. Please call create_network() first.")
            return
        
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
            with open(file_path, 'w') as f:
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

            # Convert to DataFrame and save
            results_df = pd.DataFrame(results)
            results_df.to_csv(file_path, index=False)
            
        else:
            raise ValueError(f"Unsupported format: {format}. Supported formats: gexf, gml, edgelist, csv")
        
        logger.info(f"Network saved successfully with {self.network.number_of_nodes()} nodes and {self.network.number_of_edges()} edges") 
