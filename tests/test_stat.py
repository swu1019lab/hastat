import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from itertools import combinations
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
from scipy import stats
from hastat.utils.association import HaplotypeAssociation, PAIR_COLUMNS


def frames(arrays):
    rows = [(str(i)+'_'+str(j), 'Hap'+str(i+1), x) for i, a in enumerate(arrays) for j, x in enumerate(a)]
    return (pd.DataFrame([(s,x) for s,h,x in rows], columns=['samples','trait']),
            pd.DataFrame([(s,h) for s,h,x in rows], columns=['samples','haplotypes']))


def run(arrays, method='TukeyHSD', **kwargs):
    return HaplotypeAssociation(*frames(arrays), min_hap_size=3, method=method, **kwargs)


class StatisticalValidation(unittest.TestCase):
    def close(self, observed, expected, rtol=2e-10, atol=1e-12):
        np.testing.assert_allclose(observed, expected, rtol=rtol, atol=atol)

    def test_classical_and_tukey_reference(self):
        rng = np.random.default_rng(20261028)
        examples = [[[24.5,23.5,26.4,27.1,29.9],[28.4,34.2,29.5,32.2,30.1],[26.1,28.3,24.3,26.2,27.8]]]
        examples += [[rng.normal(i*.5, 1+i*.2, n) for i,n in enumerate(ns)] for ns in [(5,5,5),(3,7,21),(156,60,7),(4,4,8,12)]]
        for arrays in examples:
            result = run(arrays)
            overall = result.result()['omnibus'][0]
            self.assertEqual(overall['status'], 'ok')
            ref = stats.f_oneway(*arrays)
            self.close([overall['statistic'],overall['pvalue']], [ref.statistic,ref.pvalue])
            allx = np.concatenate(arrays); grand = allx.mean(); k = len(arrays); n = len(allx)
            f = sum(len(a)*(np.mean(a)-grand)**2 for a in arrays)/(k-1)/(sum(sum((np.asarray(a)-np.mean(a))**2) for a in arrays)/(n-k))
            self.close(overall['statistic'], f)
            tukey = stats.tukey_hsd(*arrays); ci = tukey.confidence_interval()
            for row,(i,j) in zip(result.result()['comparisons'], combinations(range(k),2)):
                self.close(row['p-adj'], tukey.pvalue[j,i], rtol=2e-7, atol=2e-10)
                self.close([row['meandiff'],row['lower'],row['upper']],
                           [tukey.statistic[j,i],ci.low[j,i],ci.high[j,i]], rtol=2e-7, atol=2e-8)

    def test_welch_and_bonferroni_reference(self):
        rng = np.random.default_rng(807)
        for ns in [(3,5,10),(156,60,7),(3,8,4,22)]:
            arrays = [rng.normal(i*2, 1+i*3, n) for i,n in enumerate(ns)]
            result = run(arrays, 'WelchBonferroni').result(); overall = result['omnibus'][0]
            k = len(ns); n=np.array(ns); means=np.array([x.mean() for x in arrays]); var=np.array([x.var(ddof=1) for x in arrays])
            w=n/var; weighted=sum(w*means)/sum(w); b=sum((1-w/sum(w))**2/(n-1))
            f=sum(w*(means-weighted)**2)/(k-1)/(1+2*(k-2)*b/(k*k-1)); df2=(k*k-1)/(3*b)
            self.close([overall['statistic'],overall['df1'],overall['df2'],overall['pvalue']], [f,k-1,df2,stats.f.sf(f,k-1,df2)])
            m=k*(k-1)//2
            for row,(i,j) in zip(result['comparisons'], combinations(range(k),2)):
                v1,v2=var[i]/n[i],var[j]/n[j]; se=np.sqrt(v1+v2); df=(v1+v2)**2/(v1*v1/(n[i]-1)+v2*v2/(n[j]-1))
                diff=means[j]-means[i]; t=diff/se; p=2*stats.t.sf(abs(t),df)
                width=stats.t.isf(.05/(2*m),df)*se
                self.close([row['statistic'],row['df'],row['p-raw'],row['p-adj'],row['lower'],row['upper']],
                           [t,df,p,min(1,m*p),diff-width,diff+width])
                self.assertEqual(row['reject'],row['p-adj'] < .05)

    def test_two_group_welch_f_equals_t_squared(self):
        arrays=[[1,2,3,4],[2,4,7,9,14,20]]
        r=run(arrays,'WelchBonferroni').result(); pair=r['comparisons'][0]
        self.close(r['omnibus'][0]['statistic'],pair['statistic']**2)
        self.close(r['omnibus'][0]['pvalue'],pair['p-adj'])
        self.assertEqual(r['method'],'WelchBonferroni')

    def test_filter_after_join_and_per_outcome(self):
        pheno,groups=frames([[1,2,3,4],[4,5,6,7],[7,8,9,10]])
        pheno=pheno.iloc[:-1].copy(); pheno['other']=pheno.trait
        pheno.loc[pheno.samples.str.startswith('0_'),'other']=[1,np.nan,np.inf,-np.inf]
        pheno.loc[pheno.samples.str.startswith('1_'),'trait']=[4,5,np.nan,np.nan]
        r=HaplotypeAssociation(pheno,groups,min_hap_size=3).result()
        self.assertEqual([x['groups'] for x in r['omnibus']],[2,2])
        self.assertEqual([x['n'] for x in r['omnibus']],[7,7])
        self.assertEqual({x['group1'] for x in r['comparisons'] if x['pheno']=='other'},{'Hap2'})

    def test_constant_groups_retained_or_reported(self):
        arrays=[[0,0,0],[1,2,3],[3,4,7]]
        classic=run(arrays).result(); welch=run(arrays,'WelchBonferroni').result()
        self.assertEqual(classic['omnibus'][0]['groups'],3)
        self.assertEqual(len(classic['comparisons']),3)
        self.assertEqual(welch['omnibus'][0]['status'],'zero_variance')
        self.assertEqual(welch['groups'][0]['n'],3)
        self.assertEqual(welch['comparisons'],[])
        for method in ('TukeyHSD','WelchBonferroni'):
            for arrays in ([[0,0,0],[0,0,0]],[[0,0,0],[1,1,1]]):
                self.assertEqual(run(arrays,method).result()['omnibus'][0]['status'],'zero_variance')

    def test_empty_and_single_group(self):
        for arrays in ([[1,2],[4,5]],[[1,2,3],[np.nan,np.nan,np.nan]]):
            r=run(arrays)
            self.assertEqual(r.result()['omnibus'][0]['status'],'insufficient_groups')
            self.assertEqual(list(r.comparisons.columns),PAIR_COLUMNS)
            self.assertTrue(r.comparisons.empty)
            json.dumps(r.result(),allow_nan=False)

    def test_duplicate_samples_rejected(self):
        pheno,groups=frames([[1,2,3],[4,5,6]])
        for p,g in [(pd.concat([pheno,pheno.iloc[:1]]),groups),(pheno,pd.concat([groups,groups.iloc[:1]]))]:
            with self.assertRaises(ValueError): HaplotypeAssociation(p,g,min_hap_size=3)
        groups.loc[0,'samples']='  '+groups.loc[1,'samples']+' '
        with self.assertRaises(ValueError): HaplotypeAssociation(pheno,groups,min_hap_size=3)

    def test_bad_values_and_options_rejected(self):
        p,g=frames([[1,2,3],[4,5,6]])
        p=p.astype({'trait':object});p.loc[0,'trait']='bad'
        with self.assertRaises(ValueError): HaplotypeAssociation(p,g,min_hap_size=3)
        for kwargs in [{'method':'unknown'},{'method':'AllPairTest'},{'min_hap_size':1},{'alpha':0},{'alpha':float('nan')}]:
            with self.assertRaises(ValueError): HaplotypeAssociation(*frames([[1,2,3],[4,5,6]]),**kwargs)

    def test_order_invariance(self):
        p,g=frames([[1,3,5],[2,4,7,8],[5,7,9,11]])
        for method in ('TukeyHSD','WelchBonferroni'):
            a=HaplotypeAssociation(p,g,3,method=method).comparisons
            b=HaplotypeAssociation(p.sample(frac=1,random_state=8),g.sample(frac=1,random_state=3),3,method=method).comparisons
            pd.testing.assert_frame_equal(a,b,rtol=1e-12,atol=1e-12)

    def test_precision_and_csv_roundtrip(self):
        r=run([[0,1,2,3,4],[8,9,10,11,12],[4,5,7,8,9]],'WelchBonferroni')
        tiny=r.comparisons['p-adj'].min(); self.assertGreater(tiny,0); self.assertLess(tiny,.001)
        self.assertNotEqual(tiny,round(tiny,4))
        csv=r.comparisons.to_csv(index=False,float_format='%.17g')
        parsed=pd.read_csv(io.StringIO(csv),float_precision='round_trip')
        self.close(parsed['p-adj'],r.comparisons['p-adj'],rtol=0,atol=0)

    def test_tukey_decision_near_threshold_uses_unrounded_p(self):
        base=np.array([-2,-1,0,1,2.])
        critical=stats.studentized_range.ppf(.95,2,8)*np.sqrt(base.var(ddof=1)/len(base))
        for offset,expected in [(-.00001,False),(.00001,True)]:
            pair=run([base,base+critical+offset]).result()['comparisons'][0]
            self.assertEqual(round(pair['p-adj'],4),.05)
            self.assertEqual(pair['reject'],expected)

    def test_unmatched_groups_reported(self):
        p,g=frames([[1,2,3],[4,5,6]])
        p['samples']='missing_'+p.samples
        r=HaplotypeAssociation(p,g,3).result()
        self.assertEqual([x['n'] for x in r['groups']],[0,0])
        self.assertEqual(r['omnibus'][0]['status'],'insufficient_groups')

    def test_cli_empty_and_leading_zero_ids(self):
        from hastat.func.stat import run as cli_run
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);p,g=frames([[1,2,3],[4,5,7]])
            p['samples']=['00'+str(i) for i in range(6)];g['samples']=p.samples
            p.to_csv(root/'p.csv',index=False);g.to_csv(root/'g.csv',index=False)
            args=SimpleNamespace(pheno=root/'p.csv',group=root/'g.csv',min_hap_size=3,annotate='gene',method='TukeyHSD',out=root/'out.csv')
            cli_run(args)
            self.assertEqual(len(pd.read_csv(args.out)),1)
            self.assertFalse(Path(str(args.out)+'.summary.json').exists())
            args.min_hap_size=10;cli_run(args)
            self.assertTrue(pd.read_csv(args.out).empty)
            self.assertEqual(sorted(p.name for p in root.iterdir()), ['g.csv','out.csv','p.csv'])


if __name__=='__main__': unittest.main()
