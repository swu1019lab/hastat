"""Verify allele denominators, LD identity and rendered track geometry."""
import tempfile
import unittest
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
import pandas as pd
from hastat.utils.population import ld_r_squared
from hastat.viz.gene import FancyGene
from hastat.viz.tracks import allele_frequencies, aligned_genotypes
from test_gene_coordinates import coordinate_inputs
from test_population import write_vcf


class GeneFrequencyLD(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.addCleanup(plt.close, 'all')
        self.folder = Path(tmp.name)

    def test_copy_counts_partial_missing_and_population_denominators(self):
        sites = pd.DataFrame({'chrom':['1']*3,'pos':[1,2,3],'ref':'A','alt':'G',
                              'H1':[0,1,-1],'H2':[2,-1,-1]})
        a = pd.DataFrame({'samples':['001','NA','s3','s4'],'haplotypes':['H1','H1','H1','H2']})
        np.testing.assert_allclose(allele_frequencies(sites,a)['All'],[.25,.5,np.nan],equal_nan=True)
        path = self.folder/'pop.csv'
        path.write_text('sample,pop\n001,P\nNA,Q\ns3,\ns4,Q\n')
        grouped = allele_frequencies(sites,a,path)
        np.testing.assert_allclose(grouped.P,[0,.5,np.nan],equal_nan=True)
        np.testing.assert_allclose(grouped.Q,[.5,.5,np.nan],equal_nan=True)
        self.assertIsNone(allele_frequencies(sites,None))

    def test_ld_matches_exact_keys_and_entire_vcf_population(self):
        rows=[(10,'A','G',['0/0','0/1','1/1','./1','0/0']),
              (20,'A','T',['0/0','1/1','0/1','0/1','./.']),
              (30,'G','T',['0/0']*5)]
        path=write_vcf(self.folder/'input.vcf',rows,samples=list('abcde'))
        sites=pd.DataFrame({'chrom':'chr1','pos':[30,20,10],'ref':['G','A','A'],'alt':['T','T','G']})
        data=aligned_genotypes(path,sites)
        matrix=ld_r_squared(data.to_numpy())
        self.assertEqual(data.columns.tolist(),list('abcde'))
        # Only a,b,c are jointly called at 10 and 20.
        self.assertAlmostEqual(matrix[1,2],np.corrcoef([0,1,2],[0,2,1])[0,1]**2,places=6)
        self.assertTrue(np.isnan(matrix[0]).all())
        freq=allele_frequencies(sites,None,genotypes=data)
        np.testing.assert_allclose(freq.All,[0,.5,.375])
        sites.loc[1,'alt']='C'
        with self.assertRaisesRegex(ValueError,'absent from VCF'):
            aligned_genotypes(path,sites)

    def test_triangle_geometry_and_hidden_haplotypes(self):
        config,table=coordinate_inputs(self.folder)
        calls=[(p,'A','G',['0/0','0/0','1/1','0/1']) for p in table.pos]
        write_vcf(self.folder/'test.vcf',calls,chrom='1',samples=['s1','s2','s3','extra'])
        cfg=self.folder/'gene.toml'
        with cfg.open('a') as f:
            f.write('hap_cols=["Hap1"]\n')
        base=FancyGene(config); base.plot()
        def matrix_width(g):
            return g.ax.transData.transform((.9,0))[0]-g.ax.transData.transform((.1,0))[0]
        # Use automatic spacing for the windows and LD.
        cfg.write_text('\n'.join(s for s in cfg.read_text().splitlines()
                                 if not s.startswith(('gene_y=', 'hap_y=', 'ann_y=', 'pi_y=', 'fst_y=')))+'\n')
        with cfg.open('a') as f:
            f.write('layout="auto"\nvcf_file="test.vcf.gz"\nld_show=true\nshade_gene=true\n')
        obj=FancyGene(config);obj.plot();obj.fig.canvas.draw()
        mesh=next(p for p in obj.ax.collections if isinstance(p,PolyCollection) and p.get_gid()=='ld-heatmap')
        verts=np.vstack([p.vertices for p in mesh.get_paths()])
        pixels=obj.ax.transData.transform(verts)
        self.assertAlmostEqual(np.ptp(pixels[:,1])*2,np.ptp(pixels[:,0]),places=6)
        self.assertAlmostEqual(matrix_width(base),matrix_width(obj),places=5)
        self.assertEqual(len(mesh.get_array()),len(table)*(len(table)+1)//2)
        np.testing.assert_allclose(mesh.get_array(),1,atol=1e-5)
        self.assertGreater(obj.fig.get_figheight(),base.fig.get_figheight())
        # Window x axes use physical coordinates; their entire panels sit above the gene.
        gene_y=obj.ax.transData.transform(next(t for t in obj.ax.texts if t.get_text()=='G1').xy)[1]
        windows=[axis for axis in obj.ax.child_axes if axis.lines]
        self.assertEqual(len(windows),2)
        self.assertTrue(all(axis.bbox.ymin>gene_y for axis in windows))
        for axis in windows:
            shade, = [p for p in axis.patches if p.get_gid()=='gene-span']
            self.assertEqual((shade.get_x(),shade.get_width()), (99.5,101))
        # Hiding LD must not switch AF back from the full VCF to the Hap counts.
        with cfg.open('a') as f:
            f.write('[G1]\nld_show=false\naf_show=true\nshade_gene=false\n')
        hidden=FancyGene(config);hidden.plot()
        view=hidden._prepare_haplotypes(hidden.gene_data,hidden._get_gene_config('G1'))
        self.assertIsNone(view['ld'])
        np.testing.assert_allclose(view['frequencies'].All,3/8)
        self.assertFalse(any(c.get_gid()=='ld-heatmap' for c in hidden.ax.collections))
        self.assertFalse(any(p.get_gid()=='gene-span' for axis in hidden.ax.child_axes for p in axis.patches))

    def test_ld_switch_requires_vcf_and_accepts_display_controls(self):
        from hastat.viz.tracks import draw_ld
        config,_=coordinate_inputs(self.folder)
        with (self.folder/'gene.toml').open('a') as f:
            f.write('ld_show=true\n')
        with self.assertRaisesRegex(ValueError,'vcf_file'):
            FancyGene(config).plot()
        fig,ax=plt.subplots()
        draw_ld(ax,np.eye(2),0,1,1,'viridis',False)
        self.assertEqual(ax.collections[0].cmap.name,'viridis')
        self.assertFalse(ax.child_axes)

    def test_toml_positions_are_honored_with_vcf(self):
        cfg, table = coordinate_inputs(self.folder)
        calls = [(p, 'A', 'G', ['0/0','1/1']) for p in table.pos]
        write_vcf(self.folder/'test.vcf', calls, chrom='1', samples=['s1','s2'])
        path = self.folder/'gene.toml'
        path.write_text(path.read_text()+'layout="fixed"\nvcf_file="test.vcf.gz"\nld_show=true\n'
                        'ld_y=0.35\nann_legend_y=0.69\nann_band_y=0.51\n'
                        '[G1]\naf_show=true\naf_y=0.44\npi_y=0.85\nld_y=0.39\n')
        obj = FancyGene(cfg)
        obj.plot()
        np.testing.assert_allclose(obj.fig.get_size_inches(), [cfg['plot']['width'],cfg['plot']['height']])
        self.assertEqual((obj.gene_data['start'],obj.gene_data['end']), (50,220))
        pi_axis = next(a for a in obj.ax.child_axes if a.get_ylabel().startswith(r'$\mathit{\pi}$'))
        self.assertAlmostEqual(obj.ax.transData.inverted().transform(pi_axis.bbox.p0)[1], .1+.85*.8)
        alt = next(t for t in obj.ax.texts if t.get_text() == 'ALT')
        self.assertAlmostEqual(alt.get_position()[1], .1+.44*.8)
        mesh = next(c for c in obj.ax.collections if c.get_gid() == 'ld-heatmap')
        top = max(p.vertices[:,1].max() for p in mesh.get_paths())
        self.assertAlmostEqual(top, .1+.39*.8)
        bands = [p for p in obj.ax.patches if p.get_gid() == 'annotation-band']
        self.assertTrue(bands)
        for band in bands:
            self.assertAlmostEqual(band.get_y(), .1+.51*.8)
        legend = next(a for a in obj.ax.artists if a.get_gid() == 'annotation-legend')
        self.assertAlmostEqual(obj.ax.transData.inverted().transform(legend.get_bbox_to_anchor().p0)[1], .1+.69*.8)
        # VCF alone does not select a layout; the fixed coordinates also work without it.
        path.write_text(path.read_text().replace('vcf_file="test.vcf.gz"\nld_show=true\n',''))
        obj.plot()
        pi_axis = next(a for a in obj.ax.child_axes if a.get_ylabel().startswith(r'$\mathit{\pi}$'))
        self.assertAlmostEqual(obj.ax.transData.inverted().transform(pi_axis.bbox.p0)[1], .1+.85*.8)

    def test_dense_annotations_remain_visible_when_labels_cannot_fit(self):
        cfg, table = coordinate_inputs(self.folder)
        with (self.folder/'ann.csv').open('w') as stream:
            stream.write('chrom,pos,ann,label\n')
            for index in range(16):
                stream.write(f'1,100,missense,Long consequence {index}\n')
            stream.write('1,150,synonymous,\n')
        obj=FancyGene(cfg);obj.plot()
        band=[p for p in obj.ax.patches if p.get_gid()=='annotation-band']
        labels=[t for t in obj.ax.texts if t.get_text().startswith('Long consequence')]
        self.assertEqual(len(band),2)
        self.assertEqual(len(labels),2)
        renderer=obj.fig.canvas.get_renderer()
        # Annotation bboxes include leaders; compare only the drawn text boxes.
        from matplotlib.text import Text
        self.assertFalse(Text.get_window_extent(labels[0],renderer).overlaps(
            Text.get_window_extent(labels[1],renderer)))
        legend=next(a for a in obj.ax.artists if a.get_gid()=='annotation-legend')
        bounds=legend.get_window_extent(renderer)
        self.assertTrue(all(bounds.ymin>t.get_window_extent(renderer).ymax for t in labels))
        name=next(t for t in obj.ax.texts if t.get_text()=='G1')
        self.assertAlmostEqual((bounds.x0+bounds.x1)/2,obj.ax.transData.transform(name.xy)[0],places=5)
        self.assertEqual(len({t.get_window_extent(renderer).y0 for t in legend.get_texts()}),1)
