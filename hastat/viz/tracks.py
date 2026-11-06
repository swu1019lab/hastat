"""Sample-weighted allele frequencies, annotations and LD for gene haplotype figures."""
import numpy as np
import pandas as pd
from matplotlib import colors
from matplotlib.collections import PolyCollection
from matplotlib.font_manager import FontProperties
from matplotlib.legend import Legend
from matplotlib.patches import Rectangle
from matplotlib import pyplot as plt

from hastat.dataio.table import read_sample_csv, validate_sample_ids
from hastat.dataio.vcf import read as read_vcf
from hastat.log.logger import logger
from hastat.utils.variant import VariantData

KEY = ['chrom', 'pos', 'ref', 'alt']
MISSING = '#C9CED3'
REF_COLOR = '#DCE8EC'
ALT_COLOR = '#24495D'


def aligned_genotypes(path, sites):
    """Return all diploid VCF samples, ordered by the complete displayed variant keys."""
    keys = pd.MultiIndex.from_frame(sites[KEY])
    chrom = sites.chrom.iloc[0]
    with read_vcf(path) as handle:
        frame = VariantData(handle).extract_genotypes(chrom, int(sites.pos.min()), int(sites.pos.max()))
    if not frame.index.is_unique:
        raise ValueError('Duplicate VCF variant keys in the LD interval')
    missing = keys.difference(frame.index)
    if len(missing):
        raise ValueError(f'Plotted variant absent from VCF: {missing[0]} (check chrom, pos, ref, alt)')
    return frame.loc[keys]


def allele_frequencies(sites, assignments, pop_file=None, genotypes=None):
    """ALT copies / called chromosomes for each site and population."""
    if genotypes is None and assignments is None:
        logger.warning('Frequency track requires sample assignments or a VCF; skipping it')
        return None
    if genotypes is not None:
        dosages = genotypes
        units = pd.Series(genotypes.columns, index=genotypes.columns, name='unit')
    else:
        units = assignments.set_index('samples').haplotypes.rename('unit')
        dosages = sites.loc[:, units.unique()]
    if pop_file:
        pop = read_sample_csv(pop_file)
        if len(pop.columns) != 2:
            raise ValueError('pop_file requires two columns: sample and population')
        pop.columns = ['samples', 'population']
        pop = validate_sample_ids(pop).set_index('samples').population.str.strip()
        groups = pop.reindex(units.index)
        valid = groups.notna() & groups.ne('')
        if not valid.any():
            raise ValueError('No samples with population labels match the frequency input')
        if not valid.all():
            logger.warning(f'Frequency tracks exclude {(~valid).sum()} samples without a population label')
        groups, units = groups[valid], units[valid]
    else:
        groups = pd.Series('All', index=units.index)
    weights = pd.crosstab(units, groups).reindex(dosages.columns, fill_value=0)
    data = dosages.to_numpy(dtype=float)
    called = np.isfinite(data) & (data >= 0)
    alternate = np.where(called, data, 0) @ weights.to_numpy()
    chromosomes = (2 * called) @ weights.to_numpy()
    freq = np.full(alternate.shape, np.nan)
    np.divide(alternate, chromosomes, out=freq, where=chromosomes > 0)
    return pd.DataFrame(freq, index=pd.MultiIndex.from_frame(sites[KEY]), columns=weights.columns)


def draw_frequencies(ax, freq, x, y, width, height, style):
    if style not in ('bar', 'pie'):
        raise ValueError('af_style must be bar or pie')
    from matplotlib.patches import Wedge
    dx = width / len(freq)
    for i, group in enumerate(freq):
        bottom = y + (len(freq.columns) - i - 1) * height
        for j, alt in enumerate(freq[group]):
            left = x + j*dx
            if not np.isfinite(alt):
                ax.add_patch(Rectangle((left,bottom),dx,height,fc=MISSING,ec='none'))
            elif style == 'bar':
                ax.add_patch(Rectangle((left,bottom),dx,height,fc=REF_COLOR,ec='none'))
                ax.add_patch(Rectangle((left,bottom),dx,height*alt,fc=ALT_COLOR,ec='none'))
            else:
                center = (left+dx/2, bottom+height/2)
                radius = min(dx,height)*.48
                for a,b,color in [(0,360*(1-alt),REF_COLOR),(360*(1-alt),360,ALT_COLOR)]:
                    if b>a:
                        ax.add_patch(Wedge(center,radius,a,b,fc=color,ec='none'))
        ax.text(x-.012,bottom+height/2,f'{group} AF',ha='right',va='center',fontsize=7)
    ax.text(x+width+.008,y+height*.55,'REF',color='#52717F',fontsize=7)
    ax.text(x+width+.008,y,'ALT',color=ALT_COLOR,fontsize=7)


def draw_ld(ax, matrix, x, top, width, cmap_name='Blues', show_colorbar=True):
    """Rotate the pair grid 45 degrees. Each SNP owns one equal interval on top."""
    n = len(matrix)
    i,j = np.triu_indices(n)
    d = width/n
    cx = x + (i+j+1)*d/2
    cy = top - (j-i)*d/2
    vertices = np.stack([np.column_stack((cx+ox*d/2,np.minimum(top,cy+oy*d/2)))
                         for ox,oy in [(-1,0),(0,1),(1,0),(0,-1)]],axis=1)
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(MISSING)
    mesh = PolyCollection(vertices, cmap=cmap,norm=colors.Normalize(0,1),edgecolors='none',
                          antialiaseds=False, rasterized=n>250)
    mesh.set_array(np.ma.masked_invalid(matrix[i,j]))
    mesh.set_gid('ld-heatmap')
    ax.add_collection(mesh)
    if not show_colorbar:
        return
    bar = ax.inset_axes([x+.78*width,top-.31*width,.19*width,.011*width],transform=ax.transData)
    ax.figure.colorbar(mesh,cax=bar,orientation='horizontal',ticks=[0,.5,1])
    bar.tick_params(labelsize=7,length=2,pad=2)
    bar.set_title('LD ($r^2$)', fontsize=8, pad=3)
    if np.isnan(matrix).any():
        bar.text(.5,-2.5,'Gray: undefined',ha='center',va='top',fontsize=6,transform=bar.transAxes)


def read_annotations(path):
    if not path:
        return None
    data = read_sample_csv(path)
    if not {'pos','ann'} <= set(data):
        raise ValueError('Annotations require pos and ann')
    if ('ref' in data) != ('alt' in data):
        raise ValueError('Annotation ref and alt must be supplied together')
    data['pos'] = pd.to_numeric(data.pos,errors='raise')
    if not (np.isfinite(data.pos) & data.pos.ge(1) & data.pos.eq(np.floor(data.pos))).all():
        raise ValueError('Annotation positions must be positive integers')
    data['pos'] = data.pos.astype('int64')
    data['ann'] = data.ann.str.strip()
    if 'label' in data:
        data['label'] = data.label.str.strip()
    return data[data.ann.ne('')].drop_duplicates()


def draw_annotations(ax, annotations, gene, frame, sites, grid, region_to_x, *, band_bottom=None):
    """Draw a compact indexed band and return physical label anchors for final layout."""
    x,gy,w,gh = frame
    data = annotations[annotations.pos.between(gene['start'],gene['end'])]
    if 'chrom' in data:
        data = data[data.chrom.eq(gene['chrom'])]
    if data.empty:
        return [], None
    # Use two explicit semantics: ann = category, optional label = callout.
    # Allow sparse older free-text inputs to keep their individual labels.
    if 'label' not in data:
        data = data.assign(label=data.ann if len(data)<=6 else '')
    labels = [(region_to_x(int(r.pos)),gy+gh*.5,r.label) for r in data.itertuples() if r.label]
    if sites is None:
        band_x, band_y, band_w, cell = x, gy+gh*.9, w, w*.004
    else:
        band_x,band_y,band_w = grid
        cell = band_w / len(sites)
    if band_bottom is not None:
        band_y = band_bottom
    categories = data.ann.value_counts()
    names = sorted(categories.index) if len(categories)<=6 else sorted(categories.iloc[:5].index)
    misc = 'Other annotations'
    while misc in names:
        misc += '*'
    if len(categories)>6:
        logger.warning('Grouping infrequent annotation categories under Other annotations in the legend')
        names.append(misc)
    palette = ['#517D9A','#B96355','#8D9F79','#AC86A3','#C39756','#697884']
    lookup = dict(zip(names,palette))
    h = band_w*.012
    targets, physical = {}, {}
    columns = {}
    if sites is not None:
        for index, site in enumerate(sites.itertuples()):
            key = (site.pos,site.ref,site.alt) if 'ref' in data else site.pos
            columns.setdefault(key,[]).append(band_x+index*cell)
    for row in data.itertuples():
        key = (row.pos,row.ref,row.alt) if 'ref' in data else row.pos
        name = row.ann if row.ann in lookup else misc
        locations = columns.get(key,[]) if sites is not None else [region_to_x(row.pos)-cell/2]
        if not locations:
            physical.setdefault(region_to_x(row.pos),set()).add(name)
        for left in locations:
            targets.setdefault(left,set()).add(name)
    if targets:
        ax.add_patch(Rectangle((band_x,band_y),band_w,h,fc='#F0F1F2',ec='none',zorder=3))
        for left,entries in targets.items():
            for k,name in enumerate(sorted(entries)):
                patch = Rectangle((left,band_y+k*h/len(entries)),cell,h/len(entries),
                                  fc=lookup[name],ec='none',zorder=4)
                patch.set_gid('annotation-band')
                ax.add_patch(patch)
        ax.text(band_x-.008,band_y+h/2,'Annotation',ha='right',va='center',fontsize=7)
    for px, entries in physical.items():
        for k,name in enumerate(sorted(entries)):
            ax.add_patch(Rectangle((px-w*.0015, gy+gh*.55+k*h/len(entries)),w*.003,h/len(entries),
                                   fc=lookup[name],ec='none',zorder=4))
    legend = None
    if targets or physical:
        present = set().union(*targets.values(),*physical.values())
        handles = [Rectangle((0,0),1,1,fc=lookup[name]) for name in names if name in present]
        legend = Legend(ax, handles,[name for name in names if name in present],
                           loc='lower center',bbox_to_anchor=(x+w/2,gy+gh/2),ncols=len(handles),
                           bbox_transform=ax.transData,frameon=False,fontsize=7,
                           handlelength=.8,handletextpad=.45,columnspacing=1.2,borderpad=0,borderaxespad=0)
        ax.add_artist(legend)
        legend.set_gid('annotation-legend')
        legend.set_clip_on(False)
    return list(dict.fromkeys(labels)), legend


def place_labels(ax, anchors, size, label_y=None):
    """Place horizontal callouts in at most two rows, measured in physical points."""
    renderer = ax.figure.canvas.get_renderer()
    occupied = [[],[]]
    texts = []
    font = FontProperties(size=size)
    skipped = 0
    for x,y,label in sorted(anchors):
        px = ax.transData.transform((x,y))[0]
        width = renderer.get_text_width_height_descent(label,font,False)[0]+6
        interval = (px-width/2,px+width/2)
        for lane, slots in enumerate(occupied):
            if all(interval[1]<a or interval[0]>b for a,b in slots):
                slots.append(interval)
                offset = 18 if label_y is None else (
                    ax.transData.transform((x,label_y))[1]-ax.transData.transform((x,y))[1])*72/ax.figure.dpi
                text = ax.annotate(label,xy=(x,y),xytext=(0,offset+lane*(size+5)),
                            textcoords='offset points',ha='center',va='bottom',fontsize=size,
                            arrowprops=dict(arrowstyle='-',lw=.65,color='#667680'),zorder=6)
                texts.append(text)
                break
        else:
            skipped += 1
    if skipped:
        logger.warning(f'{skipped} crowded annotation labels omitted; see the annotation input for full text')
    return texts


def position_legend(ax, legend, anchor, texts):
    """Center the single legend row over the gene and clear its rendered callouts."""
    renderer = ax.figure.canvas.get_renderer()
    x_px, top = ax.transData.transform(anchor)
    for text in texts:
        if text is not None:
            top = max(top, text.get_window_extent(renderer).ymax)
    anchor = ax.transData.inverted().transform((x_px, top+4*ax.figure.dpi/72))
    legend.set_bbox_to_anchor(anchor, transform=ax.transData)
