# -*- coding: utf-8 -*-
# @Time    : 2024/1/12 17:32
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : gene.py

import numpy as np


class GeneData(object):
    """
    A class for processing results from gffutils.FeatureDB
    """

    def __init__(self, gene_id, gff_in) -> None:
        """
        :param gene_id: gene id
        :param gff_in: gffutils.FeatureDB object
        """

        self.id = gene_id
        self.gff_in = gff_in

        gene = gff_in[gene_id]
        self.chrom_name = gene.seqid
        self.chrom_start = gene.start
        self.chrom_end = gene.end
        self.id = gene.id
        self.strand = gene.strand

    def get_locus(self, upstream=0, downstream=0):
        """
        Get the locus of the gene

        :param upstream: the length of upstream with respect to the gene, default 0 bp
        :param downstream: the length of downstream with respect to the gene, default 0 bp
        :return: a tuple containing the locus of the gene
        """
        if self.strand == "+":
            start = self.chrom_start - upstream if upstream < self.chrom_start else 0
            end = self.chrom_end + downstream
        else:
            start = self.chrom_start - downstream if downstream < self.chrom_start else 0
            end = self.chrom_end + upstream
        return start, end

    def get_feature_type(self):
        """
        Get all feature types of the gff file
        :return: a list of feature type
        """
        return list(self.gff_in.featuretypes())

    def get_feature_data(self, feature_type='CDS'):
        print("Retrieving {} data of {} from FeatureDB".format(feature_type, self.id))
        return np.array([[f.start, f.end] for f in self.gff_in.region(
            seqid=self.chrom_name,
            start=self.chrom_start,
            end=self.chrom_end,
            featuretype=feature_type
        )])

    def to_gff(self, file_path=None):
        """
        Write the gene data to gff file

        :param file_path: the path of gff file
        :return: None
        """
        if file_path is None:
            file_path = "{}.gff".format(self.id)

        with open(file_path, 'w') as fn:
            for feature in self.gff_in.region(
                    seqid=self.chrom_name,
                    start=self.chrom_start,
                    end=self.chrom_end):
                fn.write(str(feature) + "\n")
