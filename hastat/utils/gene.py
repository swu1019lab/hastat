# -*- coding: utf-8 -*-
# @Time    : 2024/1/12 17:32
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : gene.py

import numpy as np
import pandas as pd
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

    def hap_table(self, chrom_name, chrom_start, chrom_end):
        """
        Get haplotypes table of a gene region

        :param chrom_name: chromosome name
        :param chrom_start: chromosome start position
        :param chrom_end: chromosome end position
        :return: a DataFrame contained haplotypes table
        """
        logger.info("Initializing the GeneVariant object")

        df = self.get_geno_data(chrom_name, chrom_start, chrom_end)
        # get unique haplotypes and their counts
        hap, counts = np.unique(df.T.values, axis=0, return_counts=True)

        # haplotypes table, including name, genotype code and size
        hap_table = pd.DataFrame(hap, columns=df.index)
        # hap_table['size'] = counts
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

    def hap_groups(self, chrom_name, chrom_start, chrom_end):
        """
        Get haplotypes group of a gene region

        :param chrom_name: chromosome name
        :param chrom_start: chromosome start position
        :param chrom_end: chromosome end position
        :return: a DataFrame contained haplotypes group
        """
        logger.info("Initializing the GeneVariant object")

        df = self.get_geno_data(chrom_name, chrom_start, chrom_end)
        # get unique haplotypes and their counts
        _, indices = np.unique(df.T.values, axis=0, return_inverse=True)

        # haplotypes group for each sample
        # add haplotype labels, e.g. Hap1, Hap2, ..., HapN
        hap_groups = pd.Series(indices, index=df.columns)
        hap_groups = "Hap" + (hap_groups + 1).astype(str)
        hap_groups.name = 'haplotypes'
        return hap_groups.reset_index()

    def calc_hap_freq(self, chrom_name, chrom_start, chrom_end, sample_groups=None):
        """
        Calculate haplotypes frequency of a gene region

        :param chrom_name: chromosome name
        :param chrom_start: chromosome start position
        :param chrom_end: chromosome end position
        :param sample_groups: A sample group Dataframe with two columns only: samples and groups
        :return: a DataFrame contained haplotypes frequency
        """
        hap_groups = self.hap_groups(chrom_name, chrom_start, chrom_end)
        if sample_groups is not None:
            return hap_groups.merge(sample_groups, how='inner', on='samples').groupby(
                ['haplotypes', 'groups']).count().unstack(fill_value=0)
        else:
            return hap_groups.groupby('haplotypes').count()

    def get_geno_data(self, chrom_name, chrom_start, chrom_end, code='code1'):
        """
        Get genotype data within a gene region

        :param chrom_name: chromosome name
        :param chrom_start: chromosome start position
        :param chrom_end: chromosome end position
        :param code: genotype code, code1 or code2, default code1, code1 is presented as 0/1/2,
                     code2 is presented as [0, 0], [0, 1] and [1, 1]
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
        return pd.DataFrame.from_dict(allele_data_set[code], orient='tight')

    def get_pi_data(self, chrom_name, chrom_start, chrom_end, size=1, step=1, sample_list=None):
        """
        Get pi data of target samples using sliding window method

        :param chrom_name: chromosome name
        :param chrom_start: chromosome start position
        :param chrom_end: chromosome end position
        :param size: window size
        :param step: step size
        :param sample_list: a list contained sample names
        :return: a numpy array contained pi data
        """
        logger.info(f"Calculating Pi data of {len(sample_list)} samples using sliding window method")
        geno_data = self.get_geno_data(chrom_name, chrom_start, chrom_end)
        if sample_list is not None:
            data = geno_data.loc[:, sample_list].to_numpy()
        else:
            data = geno_data.to_numpy()
        h = allel.HaplotypeArray(data)
        ac = h.count_alleles()
        pos = geno_data.index.get_level_values('pos').to_list()
        pi, windows, n_bases, counts = allel.windowed_diversity(pos, ac, size, step=step)
        return np.column_stack([windows, n_bases, counts, pi])

    def get_fst_data(self, chrom_name, chrom_start, chrom_end, sample_list1: list, sample_list2: list, size=1,
                     step=1):
        """
        Get Fst data of target samples using sliding window method

        :param chrom_name: chromosome name
        :param chrom_start: chromosome start position
        :param chrom_end: chromosome end position
        :param sample_list1: a list contained sample names
        :param sample_list2: a list contained sample names
        :param size: sliding window size
        :param step: step size
        :return: a numpy array contained Fst data
        """
        logger.info(
            f"Calculating Fst data between {len(sample_list1)} and {len(sample_list2)} samples "
            f"using sliding window method")
        geno_data = self.get_geno_data(chrom_name, chrom_start, chrom_end, 'code2')
        # Returns -1 for unmatched values
        idx = geno_data.columns
        index_group1 = idx.get_indexer(idx.where(idx.isin(sample_list1)).dropna())
        index_group2 = idx.get_indexer(idx.where(idx.isin(sample_list2)).dropna())
        sub_pops = [index_group1, index_group2]
        pos = geno_data.index.get_level_values('pos').to_list()
        g = allel.GenotypeArray(geno_data.to_numpy().tolist())
        _ = np.seterr(divide='ignore')
        fst, windows, counts = allel.windowed_weir_cockerham_fst(pos, g, sub_pops, size, step=step)
        return np.column_stack([windows, counts, fst])

    def get_ld_data(self, chrom_name, chrom_start, chrom_end):
        """
        Get LD data of all samples

        :return: a numpy array contained LD data
        """
        logger.info("Calculating LD data in target region {chrom_name}:{chrom_start}-{chrom_end}")
        geno_data = self.get_geno_data(chrom_name, chrom_start, chrom_end, 'code1')
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
