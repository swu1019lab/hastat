# -*- coding: utf-8 -*-
# @Time    : 2024/1/12 17:30
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : dataset.py

import os
import gffutils
import numpy as np
import pandas as pd
from collections import defaultdict
from pysam import VariantFile
from .geno import GenoData
from .gene import GeneData


class DataSet(object):
    """
    A class for loading data from VCF, GFF and CSV files
    """

    def __init__(self, vcf_file, gff_file, phe_file) -> None:
        """
        :param vcf_file: a vcf file contained genotype data
        :param gff_file: a gff file contained gene annotation data
        :param phe_file: a csv file contained phenotype data
        """
        self.vcf_in = None
        self.gff_in = None
        self.phe_in = None
        self.pheno_name = None
        self.samples_num = 0
        self.genes_num = 0
        self.pheno_num = 0
        self.set_gff(gff_file=gff_file)
        self.set_vcf(vcf_file=vcf_file)
        self.set_phe(phe_file=phe_file)

    def __str__(self) -> str:
        pars = self.samples_num, self.genes_num, self.pheno_num
        fmt = "DataSet(samples number=%d, gene number=%d, phenotypes=%d)"
        return fmt % pars

    def set_phe(self, phe_file=None):
        """Load a csv file contained phenotype data of samples

        :param phe_file: a csv file contained phenotype data of samples
        """
        if phe_file is not None:
            print("Loading the phenotype data...")
            # two columns at least, samples (fixed) and other traits name
            self.phe_in = pd.read_csv(phe_file)
            print("The number of input samples and traits from phenotype data:",
                  self.phe_in.shape[0], self.phe_in.shape[1] - 1)
            self.pheno_num = self.phe_in.shape[1] - 1
            self.pheno_name = self.phe_in.columns[1:]

    def set_gff(self, gff_file=None):
        """Load a gff file

        :param gff_file: a gff file contained gene annotation data
        """
        # Create FeatureDB object
        print("Loading the annotation data...")
        if os.path.exists(gff_file + '.sqlite3'):
            self.gff_in = gffutils.FeatureDB(gff_file + '.sqlite3')
        else:
            self.gff_in = gffutils.create_db(gff_file, gff_file + '.sqlite3')
        self.genes_num = self.gff_in.count_features_of_type(featuretype='mRNA')

    def set_vcf(self, vcf_file):
        """Load a vcf file

        :param vcf_file: a vcf file contained genotype data
        """
        print("Loading the genotype data...")
        self.vcf_in = VariantFile(vcf_file, threads=10)
        self.samples_num = len(list(self.vcf_in.header.samples))

    def get_all_samples(self):
        """Return all samples name from VCF file

        :return: a list contained all samples name
        """
        return list(self.vcf_in.header.samples)

    def get_all_genes(self):
        """Return all genes name

        :return: a generator contained all genes name
        """
        return self.gff_in.all_features(featuretype='mRNA')

    def get_gene(self, gene_id=None):
        """
        Get the gene object

        :param gene_id: gene id
        :return: GeneData object
        """
        if gene_id is not None:
            return GeneData(gene_id, self.gff_in)

    def get_feature_data(self, chrom_name, chrom_start, chrom_end, feature_type='mRNA'):
        """
        Get the feature data

        :param chrom_name: chromosome name
        :param chrom_start: chromosome start position
        :param chrom_end: chromosome end position
        :param feature_type: feature type, default is mRNA
        :return: a numpy array contained feature loci
        """
        print("Retrieving {} data of {}:{}-{} from database".format(feature_type, chrom_name, chrom_start, chrom_end))
        return np.array([[f.start, f.end] for f in self.gff_in.region(
            seqid=chrom_name,
            start=chrom_start,
            end=chrom_end,
            featuretype=feature_type
        )])

    def get_geno(self, chrom_name: str, chrom_start: int, chrom_end: int):
        """
        Get the genotype data of target region

        :param chrom_name: chromosome name
        :param chrom_start: chromosome start position
        :param chrom_end: chromosome end position
        :return: GenoData object
        """
        allele_data_set = {'code1': defaultdict(list), 'code2': defaultdict(list)}
        allele_data_set['code1']['index_names'] = ['chrom', 'pos', 'ref', 'alt']
        allele_data_set['code1']['column_names'] = ['samples']
        allele_data_set['code1']['columns'] = self.get_all_samples()
        allele_data_set['code2']['index_names'] = ['chrom', 'pos', 'ref', 'alt']
        allele_data_set['code2']['column_names'] = ['samples']
        allele_data_set['code2']['columns'] = self.get_all_samples()
        # Access snp loci
        # 0-based, half-open
        try:
            _ = next(self.vcf_in.fetch(chrom_name, chrom_start, chrom_end))
        except StopIteration:
            print("Warning: No SNP loci existed in target region!")
            return None
        except ValueError as e:
            print("Warning: {}".format(e))
            print("\tAn error occurred in {}:{}-{}".format(chrom_name, chrom_start, chrom_end))
            return None
        else:
            for rec in self.vcf_in.fetch(chrom_name, chrom_start, chrom_end):
                # print(rec.chrom, rec.pos, rec.alleles)
                # only keep bi-allelic SNP loci
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
        return GenoData(allele_data_set)

    def get_gene_geno(self, gene_id=None, upstream=0, downstream=0):
        """
        Get the genotype data of gene

        :param gene_id:gene id
        :param upstream: upstream length of gene, default is 0
        :param downstream: downstream length of gene, default is 0
        :return: GenoData object
        """
        assert gene_id is not None, "gene_id is a required input parameter"
        # access a gene
        gene = self.gff_in[gene_id]
        # get the locus of gene based on strand
        start, end = self.get_gene(gene_id).get_locus(upstream, downstream)
        print("Getting the gene genotype of {} from {}:{}-{}".format(gene.id, gene.seqid, start, end))
        gene_geno = self.get_geno(gene.seqid, start, end)
        if gene_geno:
            return gene_geno

    def get_all_genes_geno(self):
        """
        Get the genotype data of all genes

        :return: GenoData object generator
        """
        for gene in self.gff_in.all_features(featuretype='mRNA'):
            yield self.get_geno(gene.seqid, gene.start, gene.end)

    def get_promoter_geno(self, gene_id=None, promoter=2000, gene_geno=False):
        """
        Get the promoter genotype data of gene

        :param gene_id: gene id
        :param promoter: promoter length of gene, default is 2000
        :param gene_geno: whether to return gene genotype data, default is False
        :return: GenoData object
        """
        assert gene_id is not None, "gene_id is a required input parameter"
        # access a gene
        gene = self.gff_in[gene_id]
        chrom_name = gene.seqid
        if gene_geno:
            chrom_start, chrom_end = (gene.start - 1 - promoter, gene.end) if gene.strand == "+" else (
                gene.start - 1, gene.end + promoter)
            if chrom_start < 0:
                chrom_start = 0
        else:
            chrom_start, chrom_end = (gene.start - 1 - promoter, gene.start) if gene.strand == "+" else (
                gene.end, gene.end + promoter)
            if chrom_start < 0:
                chrom_start = 0
        print("Getting the promoter genotype of {} from {}:{}-{}".format(gene.id, chrom_name, chrom_start, chrom_end))
        promoter_geno = self.get_geno(chrom_name, chrom_start, chrom_end)
        if promoter_geno:
            return promoter_geno
