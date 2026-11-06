"""Haplotype profiles, counts, clustering and population comparisons."""
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial import distance
from hastat.log.logger import logger
from hastat.utils.variant import normalize_loci


def build_haplotype_tables(data, prefix='Hap'):
    """Build both outputs from the same complete, canonically ordered profiles."""
    if data.index.has_duplicates or not data.columns.is_unique:
        raise ValueError('Haplotype input must have unique variant records and sample names')
    data = data.sort_index()
    if len(data) == 0:
        return (pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt']),
                pd.DataFrame(columns=['samples', 'haplotypes']))
    if data.shape[1] == 0 or not data.isin([0, 1, 2]).all().all():
        raise ValueError('Haplotype grouping requires samples with complete 0/1/2 states')
    hap, inverse, counts = np.unique(data.T.to_numpy(dtype=np.int8), axis=0,
                                     return_inverse=True, return_counts=True)
    order = sorted(range(len(hap)), key=lambda i: (int(counts[i]), tuple(hap[i])), reverse=True)
    ranks = np.empty(len(order), dtype=int)
    ranks[order] = np.arange(len(order))
    labels = [f'{prefix}{i + 1}' for i in range(len(order))]
    table = pd.DataFrame(hap[order].T, index=data.index, columns=labels).reset_index()
    groups = pd.DataFrame({'samples': data.columns.tolist(),
                           'haplotypes': [labels[ranks[i]] for i in inverse]})
    return table, groups


def count_haplotypes(assignments, sample_groups, available_samples):
    """Count assigned individuals, optionally restricted to a validated group file."""
    order = sorted(assignments.haplotypes.unique(), key=lambda hap: int(hap[3:]))
    if sample_groups is None:
        counts = assignments.groupby('haplotypes').size().reindex(order, fill_value=0)
        counts.index.name = 'haplotypes'
        return counts.astype('int64').rename('samples').to_frame()
    if (not isinstance(sample_groups, pd.DataFrame) or len(sample_groups.columns) != 2
            or set(sample_groups.columns) != {'samples', 'groups'}):
        raise ValueError('Frequency groups must contain exactly samples and groups columns')
    metadata = sample_groups.loc[:, ['samples', 'groups']].copy()
    if metadata.empty or any(not isinstance(value, str) or not value.strip()
                             for value in metadata.to_numpy().flat):
        raise ValueError('Frequency groups require non-empty string sample and group names')
    if metadata.samples.duplicated().any():
        raise ValueError('Frequency group file must assign each sample exactly once')
    unknown = metadata.loc[~metadata.samples.isin(available_samples), 'samples'].tolist()
    if unknown:
        raise ValueError(f"Samples not found in VCF: {', '.join(unknown)}")
    labels = sorted(metadata.groups.unique())
    if 'haplotypes' in labels:
        raise ValueError("Population label 'haplotypes' is reserved for the frequency row identifier")
    if len(metadata) < len(available_samples):
        logger.info(f'Frequency counts use the {len(metadata)} listed samples; '
                    'haplotype definitions and ranks still use all VCF samples')
    joined = assignments.merge(metadata, on='samples', how='inner', validate='one_to_one')
    if joined.empty:
        return pd.DataFrame(index=pd.Index([], name='haplotypes'), columns=labels, dtype='int64')
    counts = joined.groupby(['haplotypes', 'groups']).size().unstack(fill_value=0)
    observed = [hap for hap in order if hap in counts.index]
    counts = counts.reindex(index=observed, columns=labels, fill_value=0).astype('int64')
    counts.index.name, counts.columns.name = 'haplotypes', None
    return counts


def extract_haplotypes(variants, loci, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
    """Return a profile table and sample groups for 1-based inclusive intervals."""
    loci = normalize_loci(loci)
    options = dict(het_threshold=het_threshold, maf_threshold=maf_threshold,
                   max_missing_threshold=max_missing_threshold)
    data = (variants.extract_genotypes(*loci[0], **options) if len(loci) == 1
            else variants.extract_multi_loci_genotypes(loci, **options))
    if data.index.has_duplicates:
        raise ValueError("Duplicate variant records in haplotype input")
    complete = data.isin([0, 1, 2]).all(axis=1)
    excluded = int((~complete).sum())
    if excluded:
        logger.warning(f"Haplotype grouping excludes {excluded} of {len(data)} records with missing calls; "
                       "all samples use the same complete loci")
    return build_haplotype_tables(data.loc[complete], 'Hap')


def extract_haplotype_table(variants, chrom_name, chrom_start, chrom_end, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
    """Return profiles on complete loci, retaining the original table columns."""
    return extract_haplotypes(variants, [(chrom_name, chrom_start, chrom_end)], het_threshold, maf_threshold, max_missing_threshold)[0]


def extract_multi_loci_haplotype_table(variants, loci, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
    """Return profiles across the union of requested intervals."""
    return extract_haplotypes(variants, loci, het_threshold, maf_threshold, max_missing_threshold)[0]


def extract_haplotype_groups(variants, chrom_name, chrom_start, chrom_end, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
    """Assign samples using the same profiles/ranks as extract_haplotype_table."""
    return extract_haplotypes(variants, [(chrom_name, chrom_start, chrom_end)], het_threshold, maf_threshold, max_missing_threshold)[1]


def extract_multi_loci_haplotype_groups(variants, loci, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
    """Assign samples across the union of requested intervals."""
    return extract_haplotypes(variants, loci, het_threshold, maf_threshold, max_missing_threshold)[1]


def calculate_haplotype_counts(variants, chrom_name, chrom_start, chrom_end, sample_groups=None, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
    """Return integer sample counts per profile; group-file subsets affect counts only."""
    groups = extract_haplotype_groups(variants, chrom_name, chrom_start, chrom_end, het_threshold, maf_threshold, max_missing_threshold)
    return count_haplotypes(groups, sample_groups, list(variants.vcf_in.header.samples))


def calculate_multi_loci_haplotype_counts(variants, loci, sample_groups=None, het_threshold=None, maf_threshold=0.0, max_missing_threshold=1.0):
    """Count multi-region profile assignments with the same rule as single-region freq."""
    groups = extract_multi_loci_haplotype_groups(variants, loci, het_threshold, maf_threshold, max_missing_threshold)
    return count_haplotypes(groups, sample_groups, list(variants.vcf_in.header.samples))


def cluster_haplotypes(table, haplotypes, method='average'):
    """Return ordered IDs and linkage on common complete loci (Hamming distance).

    Stable ID order resolves input-order dependence in distance ties. Each profile
    has equal weight; sample counts affect labels only. This is a display tree,
    not a phylogenetic inference or statistical significance test.
    """
    allowed = ('single', 'complete', 'average', 'weighted')
    if method not in allowed:
        raise ValueError('hap_cluster_linkage must be single, complete, average or weighted '
                         '(Hamming distance); ward, centroid and median require Euclidean distances')
    names = sorted(haplotypes)
    if len(set(names)) != len(names):
        raise ValueError('Clustering requires unique haplotype IDs')
    values = table[names].to_numpy(dtype=float)
    if not np.isin(values, [-1, 0, 1, 2]).all():
        raise ValueError('Clustering genotypes must be -1 (missing), 0, 1 or 2')
    if len(names) < 2:
        return names, None
    complete = (values != -1).all(axis=1)
    excluded = int((~complete).sum())
    if excluded:
        logger.warning(f'Haplotype clustering excludes {excluded} loci with missing calls; '
                       'all pairs use the same complete loci')
    if not complete.any():
        raise ValueError('No complete loci shared by displayed haplotypes for clustering')
    distances = distance.pdist(values[complete].T, metric='hamming')
    linkage = hierarchy.linkage(distances, method=method)
    leaves = hierarchy.leaves_list(linkage)
    return [names[i] for i in leaves], linkage


class HaplotypeComparison(object):
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
            data = (gv.extract_genotypes(*loci[0]) if len(loci) == 1
                    else gv.extract_multi_loci_genotypes(loci)).copy()
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
        return build_haplotype_tables(data, prefix)

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
