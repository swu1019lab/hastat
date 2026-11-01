"""One-way mean comparisons; see STAT_METHODS.md for formulas and validation."""
from itertools import combinations
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.oneway import anova_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd

from hastat import __version__ as STAT_VERSION
from hastat.utils.sample_data import validate_sample_ids, numeric_trait
METHODS = ('TukeyHSD', 'WelchBonferroni')
PAIR_COLUMNS = ['pheno', 'annotate', 'method', 'group1', 'group2', 'count1', 'count2',
                'mean1', 'mean2', 'std1', 'std2', 'meandiff', 'lower', 'upper',
                'p-raw', 'p-adj', 'reject', 'statistic', 'df', 'anova', 'family_size']


def _validated_frame(frame, name):
    return validate_sample_ids(frame)


class HapAnovaTest:
    def __init__(self, pheno, groups, min_hap_size=10, annotate='gene', method='TukeyHSD', alpha=0.05):
        if method not in METHODS:
            raise ValueError(f'Unknown comparison method {method!r}; choose TukeyHSD or WelchBonferroni')
        if not isinstance(min_hap_size, (int, np.integer)) or isinstance(min_hap_size, bool) or min_hap_size < 2:
            raise ValueError('min_hap_size must be an integer >= 2')
        if not 0 < alpha < 1:
            raise ValueError('alpha must lie strictly between zero and one')
        pheno = _validated_frame(pheno, 'Phenotypes')
        groups = _validated_frame(groups, 'Groups')
        if set(groups.columns) != {'samples', 'haplotypes'} or 'haplotypes' in pheno or len(pheno.columns) < 2:
            raise ValueError('Expected phenotype columns and samples/haplotypes group columns')
        if groups['haplotypes'].isna().any():
            raise ValueError('Missing haplotype assignment')
        groups['haplotypes'] = groups['haplotypes'].astype(str).str.strip()
        if (groups['haplotypes'] == '').any():
            raise ValueError('Empty haplotype assignment')
        traits = [c for c in pheno if c != 'samples']
        for trait in traits:
            pheno[trait] = numeric_trait(pheno[trait])
        merged = pheno.merge(groups, how='inner', on='samples', validate='one_to_one')
        self.group_labels = sorted(groups['haplotypes'].unique())
        self.phap_data = merged.set_index(['samples', 'haplotypes'])[traits]
        self._stat_des = self.phap_data.groupby(level='haplotypes').describe()
        self._stat_anova = {'comparisons': [], 'omnibus': [], 'groups': []}
        self.alpha, self.min_hap_size = float(alpha), int(min_hap_size)
        self.method = method
        self.annotate = annotate
        for trait in traits:
            self.perform_anova_and_comparisons(trait)

    def perform_anova_and_comparisons(self, phe):
        data = self.phap_data[phe].dropna().reset_index()
        all_groups = self.group_labels
        arrays, labels = [], []
        for label in all_groups:
            values = data.loc[data.haplotypes == label, phe].to_numpy(dtype=float)
            n = len(values)
            eligible = n >= self.min_hap_size
            self._stat_anova['groups'].append(dict(pheno=phe, group=label, n=n,
                mean=float(values.mean()) if n else None,
                std=float(values.std(ddof=1)) if n > 1 else None,
                status='included' if eligible else 'insufficient_samples'))
            if eligible:
                labels.append(label)
                arrays.append(values)
        k, n = len(arrays), sum(map(len, arrays))
        equal = self.method == 'TukeyHSD'
        overall = dict(pheno=phe, annotate=self.annotate, method=self.method,
            test='One-way ANOVA' if equal else 'Welch ANOVA',
            status='insufficient_groups', statistic=None, df1=None, df2=None, pvalue=None,
            groups=k, n=n, alpha=self.alpha, min_hap_size=self.min_hap_size,
            family_size=k*(k-1)//2, variance='equal' if equal else 'unequal')
        self._stat_anova['omnibus'].append(overall)
        if k < 2:
            return
        variances = np.array([a.var(ddof=1) for a in arrays])
        if (not equal and np.any(variances == 0)) or np.all(variances == 0):
            overall['status'] = 'zero_variance'
            return
        result = anova_oneway(arrays, use_var='equal' if equal else 'unequal', welch_correction=True)
        if not all(np.isfinite(v) for v in (result.statistic, result.pvalue, *result.df)):
            overall['status'] = 'not_estimable'
            return
        overall.update(status='ok', statistic=float(result.statistic), pvalue=float(result.pvalue),
                       df1=float(result.df[0]), df2=float(result.df[1]))
        m = overall['family_size']
        tukey = None
        if equal:
            tukey = pairwise_tukeyhsd(np.concatenate(arrays), np.repeat(labels, list(map(len, arrays))), alpha=self.alpha)
            if list(tukey.groupsunique) != labels:
                raise RuntimeError('Unexpected Tukey group ordering')
        rows = []
        for index, (i, j) in enumerate(combinations(range(k), 2)):
            a, b = arrays[i], arrays[j]
            row = dict(pheno=phe, annotate=self.annotate, method=self.method,
                group1=labels[i], group2=labels[j], count1=len(a), count2=len(b),
                mean1=float(a.mean()), mean2=float(b.mean()),
                std1=float(a.std(ddof=1)), std2=float(b.std(ddof=1)),
                meandiff=float(b.mean()-a.mean()), anova=overall['pvalue'], family_size=m,
                **{'p-raw': None, 'statistic': None, 'df': None})
            if equal:
                row.update(lower=float(tukey.confint[index, 0]), upper=float(tukey.confint[index, 1]),
                           reject=bool(tukey.reject[index]), **{'p-adj': float(tukey.pvalues[index])})
            else:
                test = stats.ttest_ind(b, a, equal_var=False)
                va, vb = variances[i]/len(a), variances[j]/len(b)
                df = (va+vb)**2/(va**2/(len(a)-1)+vb**2/(len(b)-1))
                half_width = stats.t.ppf(1-self.alpha/(2*m), df)*np.sqrt(va+vb)
                adjusted = min(1.0, float(test.pvalue)*m)
                row.update(lower=float(row['meandiff']-half_width), upper=float(row['meandiff']+half_width),
                    statistic=float(test.statistic), df=float(df), reject=adjusted < self.alpha,
                    **{'p-raw': float(test.pvalue), 'p-adj': adjusted})
            rows.append(row)
        self._stat_anova['comparisons'].append(pd.DataFrame(rows, columns=PAIR_COLUMNS))

    @property
    def basic_data(self):
        return self._stat_des

    @property
    def anova_data(self):
        return self._stat_anova

    @property
    def comparisons(self):
        frames = self._stat_anova['comparisons']
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=PAIR_COLUMNS)

    def result(self):
        return dict(version=STAT_VERSION, method=self.method, alpha=self.alpha,
                    omnibus=self._stat_anova['omnibus'], groups=self._stat_anova['groups'],
                    comparisons=self.comparisons.to_dict(orient='records'))
