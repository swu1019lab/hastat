# -*- coding: utf-8 -*-
# @Time    : 2024/9/30 18:14
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : gff.py

import os
import gffutils
import numpy as np
from hastat.log.logger import logger


def read(file_path: str) -> gffutils.FeatureDB:
    """
    Load a gff file contained gene annotation data of samples

    :param file_path: the path of the gff file
    :return: a FeatureDB object
    """
    logger.info("Reading the gff file: %s", file_path)
    try:
        if os.path.exists(file_path + '.sqlite3'):
            gff_in = gffutils.FeatureDB(file_path + '.sqlite3')
        else:
            # Use merge_strategy='create_unique' to handle duplicate IDs (e.g. exons/CDS with same ID)
            # by appending a unique suffix (e.g. _1, _2)
            gff_in = gffutils.create_db(file_path, dbfn=file_path + '.sqlite3', merge_strategy='create_unique')
        for feature in gff_in.featuretypes():
            logger.info("The number of %s in the gff file: %d", feature, gff_in.count_features_of_type(feature))
        return gff_in
    except Exception as e:
        logger.error("Failed to read the gff file: %s", file_path)
        raise e


class GeneAnnotation(object):
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
