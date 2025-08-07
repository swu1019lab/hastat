# -*- coding: utf-8 -*-
# @Time    : 2024/9/30 18:07
# @Author  : LXD
# @Email   : lxd1997xy@163.com
# @File    : view.py

import time
import pandas as pd
from itertools import combinations
from hastat.log.logger import logger

# Set the display option to show all rows
pd.set_option('display.max_rows', None)

# Set the display option to show all columns
pd.set_option('display.max_columns', None)


def run(args):
    logger.info('Start to run view function')

    from hastat.utils import gene
    from hastat.dataio import gff, vcf

    start_time = time.time()
    
    # Handle compare type specially
    if args.type == 'compare':
        if len(args.vcf) < 2:
            logger.error("Compare analysis requires at least 2 VCF files!")
            raise ValueError("Compare analysis requires at least 2 VCF files!")
        
        if args.pop_names and len(args.pop_names) != len(args.vcf):
            logger.error("Number of population names must match number of VCF files!")
            raise ValueError("Number of population names must match number of VCF files!")
        
        # Use default population names if not provided
        pop_names = args.pop_names if args.pop_names else [f"Pop{i+1}" for i in range(len(args.vcf))]
        
        # Create GeneVariant objects for each VCF
        gv_list = [gene.GeneVariant(vcf.read(vcf_file)) for vcf_file in args.vcf]
    else:
        # For single VCF file analysis, use the first VCF file
        gv = gene.GeneVariant(vcf.read(args.vcf[0]))

    # Initialize variables
    multi_loci = False
    merged_loci = []
    
    if args.region:
        logger.info(f"Analyzing region: {args.region}")
        chrom, start, end = args.region.replace(' ', '').replace('-', ':').split(':')
        merged_loci = [(chrom, int(start), int(end))]
    elif args.gene_id:
        logger.info(f"Analyzing gene: {args.gene_id}")
        gf = gene.GeneFeature(gff.read(args.gff))
        chrom, start, end = gf.get_locus(args.gene_id, args.upstream, args.downstream)
        merged_loci = [(chrom, start, end)]
    elif args.homo:
        if not args.gff:
            logger.error("The GFF file should be provided when using --homo!!!")
            raise ValueError("The GFF file should be provided when using --homo!!!")
        logger.info(f"Analyzing homologous genes: {args.homo}")
        gf = gene.GeneFeature(gff.read(args.gff))
        gene_ids = [gene_id.strip() for gene_id in args.homo.split(',')]
        merged_loci = gf.get_merged_locus(gene_ids, args.upstream, args.downstream)
        multi_loci = True
    else:
        logger.error("The region, gene ID, or homologous genes should be provided!!!")
        raise ValueError("The region, gene ID, or homologous genes should be provided!!!")

    data = pd.DataFrame()
    if args.type == 'compare':
        # Multi-population haplotype comparison analysis: 2025-07-25
        logger.info(f"Performing multi-population comparison with {len(gv_list)} populations: {', '.join(pop_names)}")
        
        # Create a MultiPopulationAnalyzer instance
        analyzer = gene.MultiPopulationAnalyzer(gv_list, pop_names)
        
        # Perform comparison analysis
        results = analyzer.compare_populations(merged_loci if multi_loci else [merged_loci[0]])
        
        # Save results to multiple files
        if args.out:
            # Population haplotype table and groups based on core variants
            for pop_name in pop_names:
                if f'{pop_name}_hap_table' in results and not results[f'{pop_name}_hap_table'].empty:
                    results[f'{pop_name}_hap_table'].to_csv(f"{args.out}_{pop_name}.table.csv", index=False)
                    results[f'{pop_name}_hap_groups'].to_csv(f"{args.out}_{pop_name}.group.csv", index=False)
                    logger.info(f"Population {pop_name} haplotype results saved: {args.out}_{pop_name}.table.csv and {args.out}_{pop_name}.group.csv")
            
            # Merged populations haplotype table and groups (based on core variants)
            if not results['merged_hap_table'].empty:
                results['merged_hap_table'].to_csv(f"{args.out}.merged.table.csv", index=False)
                results['merged_hap_groups'].to_csv(f"{args.out}.merged.group.csv", index=False)
                logger.info(f"Merged population results saved: {args.out}.merged.table.csv and {args.out}.merged.group.csv")
            
            # Summary comparison results
            results['comparison_summary'].to_csv(f"{args.out}.comparison_summary.csv", index=False)
            results['variant_class'].to_csv(f"{args.out}.variant_class.csv", index=False)
            
            logger.info(f"Comparison analysis completed.")
        else:
            # Print summary to stdout
            print(results['comparison_summary'])
    elif args.type == 'geno':
        data = gv.get_multi_loci_geno_data(merged_loci, het_threshold=args.het) if multi_loci else gv.get_geno_data(*merged_loci[0], het_threshold=args.het)
        if data.empty:
            logger.warning(f"The {args.type} data of target gene or region is empty!!!")
            return
        if args.out:
            data.to_csv(args.out + '.' + args.type + '.csv', index=False)
        else:
            print(data)
    elif args.type == 'table':
        data = gv.multi_loci_hap_table(merged_loci, het_threshold=args.het) if multi_loci else gv.hap_table(*merged_loci[0], het_threshold=args.het)
        if data.empty:
            logger.warning(f"The {args.type} data of target gene or region is empty!!!")
            return
        if args.out:
            data.to_csv(args.out + '.' + args.type + '.csv', index=False)
        else:
            print(data)
    elif args.type == 'group':
        data = gv.multi_loci_hap_groups(merged_loci, het_threshold=args.het) if multi_loci else gv.hap_groups(*merged_loci[0], het_threshold=args.het)
        if data.empty:
            logger.warning(f"The {args.type} data of target gene or region is empty!!!")
            return
        if args.out:
            data.to_csv(args.out + '.' + args.type + '.csv', index=False)
        else:
            print(data)
    elif args.type == 'freq':
        groups = pd.read_csv(args.group, header=0, names=['samples', 'groups']) if args.group else None
        if multi_loci:
            data = gv.calc_multi_loci_hap_freq(merged_loci, groups, het_threshold=args.het)
        else:
            data = gv.calc_hap_freq(*merged_loci[0], groups, het_threshold=args.het)
        data.reset_index(inplace=True)
        if data.empty:
            logger.warning(f"The {args.type} data of target gene or region is empty!!!")
            return
        if args.out:
            data.to_csv(args.out + '.' + args.type + '.csv', index=False)
        else:
            print(data)
    elif args.type == 'pi':
        if not args.group:
            logger.error("The group file should be provided!!!")
            raise ValueError("The group file should be provided!!!")
        groups = pd.read_csv(args.group, header=0, names=['samples', 'groups'])
        data_list = []
        # list all groups to calculate Pi
        for group in groups.groups.unique():
            if multi_loci:
                pi_data = gv.get_multi_loci_pi_data(merged_loci,
                                      size=getattr(args, 'size', 1),
                                      step=getattr(args, 'step', 1),
                                      sample_list=groups.query('groups == @group').samples.tolist(),
                                      het_threshold=args.het)
            else:
                pi_data = gv.get_pi_data(*merged_loci[0],
                                      size=getattr(args, 'size', 1),
                                      step=getattr(args, 'step', 1),
                                      sample_list=groups.query('groups == @group').samples.tolist(),
                                      het_threshold=args.het)
            data_list.append(pd.DataFrame(pi_data, columns=['chrom', 'start', 'end', 'n_bases', 'counts', 'pi']).assign(groups=group))
        data = pd.concat(data_list)
        if data.empty:
            logger.warning(f"The {args.type} data of target gene or region is empty!!!")
            return
        if args.out:
            data.to_csv(args.out + '.' + args.type + '.csv', index=False)
        else:
            print(data)
    elif args.type == 'fst':
        if not args.group:
            logger.error("The group file should be provided!!!")
            raise ValueError("The group file should be provided!!!")
        groups = pd.read_csv(args.group, header=0, names=['samples', 'groups'])
        data_list = []
        # list all pairwise comparisons to calculate Fst
        for pair in combinations(groups.groups.unique(), 2):
            if multi_loci:
                fst_data = gv.get_multi_loci_fst_data(merged_loci,
                                       groups.query('groups == @pair[0]').samples.tolist(),
                                       groups.query('groups == @pair[1]').samples.tolist(),
                                       size=getattr(args, 'size', 1),
                                       step=getattr(args, 'step', 1),
                                       het_threshold=args.het)
            else:
                fst_data = gv.get_fst_data(*merged_loci[0],
                                       groups.query('groups == @pair[0]').samples.tolist(),
                                       groups.query('groups == @pair[1]').samples.tolist(),
                                       size=getattr(args, 'size', 1),
                                       step=getattr(args, 'step', 1),
                                       het_threshold=args.het)
            data_list.append(pd.DataFrame(fst_data, columns=['chrom', 'start', 'end', 'counts', 'fst']).assign(groups='_'.join(pair)))
        data = pd.concat(data_list)
        if data.empty:
            logger.warning(f"The {args.type} data of target gene or region is empty!!!")
            return
        if args.out:
            data.to_csv(args.out + '.' + args.type + '.csv', index=False)
        else:
            print(data)
    else:
        logger.error("The data type {} is not supported!!!".format(args.type))

    end_time = time.time()
    logger.info("view function runs {:.2f} seconds".format(end_time - start_time))
