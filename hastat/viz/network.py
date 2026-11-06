"""Plot one haplotype network with sample composition and count legends."""

import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.font_manager import FontProperties
from matplotlib.offsetbox import AnnotationBbox, AnchoredOffsetbox, DrawingArea, HPacker, TextArea, VPacker
from matplotlib.transforms import Bbox
import networkx as nx
import numpy as np
import pandas as pd

from hastat.dataio.network import read_network
from hastat.dataio.table import merge_samples, normalize_sample_table, read_sample_csv, validate_sample_ids
from hastat.log.logger import logger
from hastat.viz._export import save_figure


class FancyNetwork:
    def __init__(self, config):
        self.config = config
        data = config['data']
        self.network = read_network(data['file'])
        if not len(self.network):
            raise ValueError("An empty network cannot be plotted")
        for u, v, edge in self.network.edges(data=True):
            weight = float(edge.get('weight', 1))
            if u == v or not np.isfinite(weight) or weight < 0:
                raise ValueError("Network edges must join distinct nodes with finite nonnegative weights")
            edge['weight'] = weight
        self.sample_hap_data = None
        self.hap_composition = {}
        self.group_colors = {}
        self.node_sizes = dict.fromkeys(self.network, 1)
        if not data.get('sample_hap'):
            if data.get('sample_group'):
                raise ValueError("A sample group file requires sample haplotype assignments")
            return
        self.sample_hap_data = normalize_sample_table(read_sample_csv(data['sample_hap']), 'haplotypes')
        samples = self.sample_hap_data
        if data.get('sample_group'):
            meta = read_sample_csv(data['sample_group'])
            index = data.get('group_index', 1)
            if type(index) is not int or not 1 <= index < len(meta.columns):
                raise ValueError(f"group_index must be between 1 and {len(meta.columns)-1}")
            self.group_label = str(meta.columns[index])
            selected = meta.iloc[:, [0, index]].copy()
            selected.columns = ['sample', 'population']
            selected = validate_sample_ids(selected, 'sample')
            selected.population = selected.population.str.strip()
            samples = merge_samples(samples, selected)
            valid = samples.population.notna() & samples.population.ne('')
            if not valid.all():
                logger.warning("Excluded %d samples without a selected group", int((~valid).sum()))
            samples = samples.loc[valid]
            if samples.empty:
                raise ValueError("No samples with a nonempty selected group")
            matched = samples[samples.haplotypes.isin(self.network)]
            counts = pd.crosstab(matched.haplotypes, matched.population)
            self.group_colors = self._generate_colors(sorted(counts.columns), config['plot'].get('node_color'))
            for node, row in counts.iterrows():
                self.hap_composition[node] = (row[row > 0] / row.sum()).to_dict()
        counts = samples.haplotypes.value_counts()
        self.node_sizes = {n:int(counts.get(n, 0)) for n in self.network}

    def _generate_colors(self, groups, custom_colors=None):
        cmap = LinearSegmentedColormap.from_list('hastat_groups',
            ['#C5504B', '#114F8B', '#FCE988', '#90CAEE'], N=100)
        colors = dict(zip(groups, cmap(np.linspace(0, 1, len(groups)))[:, :3]))
        for group, value in zip(groups, custom_colors or []):
            colors[group] = to_rgb(value)
        return colors

    def _area_scale(self):
        scale = self.config['plot'].get('node_scale_factor')
        if scale is None:
            diameter = 52 if self.sample_hap_data is not None else 28
            return diameter**2 / max(1, max(self.node_sizes.values(), default=1))
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("node_scale_factor must be finite and positive")
        return scale

    def _size_legend_counts(self):
        if self.sample_hap_data is None:
            return []
        counts = sorted({s for s in self.node_sizes.values() if s > 0})
        return counts if len(counts) <= 3 else [counts[0], counts[len(counts)//2], counts[-1]]

    def _get_layout_positions(self, layout='kamada_kawai', seed=None, k=None):
        # Stable graph order makes the same input independent of CSV row order.
        g = nx.Graph()
        g.add_nodes_from(sorted(self.network))
        g.add_weighted_edges_from(sorted((min(u,v), max(u,v), e['weight'])
                                         for u,v,e in self.network.edges(data=True)))
        if layout == 'kamada_kawai':
            if len(g) < 2:
                return {node:np.zeros(2) for node in g}
            distances = {a:{b:(max(float(d), .5) if a != b else 0)
                             for b,d in row.items()}
                         for a,row in nx.all_pairs_dijkstra_path_length(g, weight='weight')}
            return nx.kamada_kawai_layout(g, dist=distances)
        if layout in ('spring', 'fruchterman_reingold'):
            return nx.spring_layout(g, seed=seed, k=k)
        if layout == 'random':
            return nx.random_layout(g, seed=seed)
        return {'circular':nx.circular_layout, 'shell':nx.shell_layout,
                'spectral':nx.spectral_layout}[layout](g)

    def _space_nodes(self, raw, radii, width, height):
        names = sorted(raw)
        points = np.array([raw[n] for n in names], dtype=float)
        if len(points) < 2:
            return {n:np.array([width/2,height/2]) for n in names}
        points -= points.mean(axis=0)
        _, _, axes = np.linalg.svd(points, full_matrices=False)
        points = points @ axes.T
        # Fix the otherwise arbitrary signs of the principal axes.
        for i, ascending in enumerate((True, False)):
            if (points[:,i].argmin() < points[:,i].argmax()) != ascending:
                points[:,i] *= -1
        extents = np.ptp(points, axis=0)
        padding = 2*max(radii.values()) + 18
        usable = np.maximum([width-padding, height-padding], 1)
        points *= min(usable / np.maximum(extents, 1e-9))
        for _ in range(150):
            moved = 0.
            for i,u in enumerate(names):
                for j in range(i+1,len(names)):
                    v = names[j]
                    delta = points[j]-points[i]
                    distance = np.linalg.norm(delta)
                    gap = 24 if self.network.has_edge(u,v) else 6
                    shift = (radii[u]+radii[v]+gap-distance)/2
                    if shift > 0:
                        direction = delta/distance if distance else np.array([1.,0.])
                        points[i] -= shift*direction
                        points[j] += shift*direction
                        moved = max(moved,shift)
            if moved < .001:
                break
        points -= (points.min(axis=0)+points.max(axis=0))/2
        return dict(zip(names,points + [width/2,height/2]))

    def _draw_pie_node(self, ax, pos, composition, size, scale_factor, node_alpha):
        radius = np.sqrt(size*scale_factor)/2
        if radius == 0:
            return
        box = DrawingArea(2*radius, 2*radius)
        center = (radius,radius)
        box.add_artist(patches.Circle(center,radius,fc=ax.get_facecolor(),ec='none'))
        angle = 90.
        if composition:
            for group, fraction in sorted(composition.items()):
                stop = angle+360*fraction
                box.add_artist(patches.Wedge(center,radius,angle,stop,fc=self.group_colors[group],
                                             ec='none',alpha=node_alpha))
                angle=stop
        else:
            box.add_artist(patches.Circle(center,radius,fc='C3',ec='none',alpha=node_alpha))
        box.add_artist(patches.Circle(center,radius,fill=False,ec='white',lw=.55))
        artist = AnnotationBbox(box,pos,xycoords='data',frameon=False,pad=0,annotation_clip=False)
        artist.set_zorder(4)
        ax.add_artist(artist)

    def _legend(self, ax):
        sections = []
        cfg = self.config['plot']
        def text(label, bold=False):
            return TextArea(label,textprops=dict(size=9,weight='semibold' if bold else 'normal'))
        def key(label, color):
            area = DrawingArea(9,10)
            area.add_artist(patches.Rectangle((0,1),8,8,fc=color,ec='none'))
            return HPacker(children=[area,text(label)],pad=0,sep=5,align='center')
        if self.group_colors:
            sections.append(VPacker(children=[text(getattr(self,'group_label','Group').replace('_',' ').capitalize(),True)]+[
                key(group,color) for group,color in self.group_colors.items()],pad=0,sep=3,align='left'))
        counts = self._size_legend_counts()
        self.size_legend_markers, self.size_legend_labels = [], []
        if counts:
            radii = np.sqrt(np.array(counts)*self._area_scale())/2
            largest = radii[-1]
            ys = []
            for r in radii:
                gap = 9+max(0,cfg.get('legend_label_spacing',2))
                ys.append(max(2*r,ys[-1]+gap if ys else 2*r))
            label_x = 2*largest + 8 + max(0,cfg.get('legend_handle_text_pad',1))*3
            area = DrawingArea(label_x + 6*max(len(str(n)) for n in counts),ys[-1]+6)
            for n,r,y in zip(counts,radii,ys):
                circle = patches.Circle((largest,r),r,fill=False,ec='#525D65',lw=.85)
                area.add_artist(circle)
                area.add_artist(plt.Line2D([largest+1,label_x-3],[2*r,y],color='#899198',lw=.55,ls=':'))
                label = plt.Text(label_x,y,str(n),fontsize=8,va='center')
                area.add_artist(label)
                self.size_legend_markers.append(circle)
                self.size_legend_labels.append(label)
            sections.append(VPacker(children=[text('Samples (n)',True),area],sep=3,pad=0,align='left'))
        if not sections:
            return 0.
        self.legend_box = AnchoredOffsetbox(loc='lower center',pad=0,borderpad=0,
            child=HPacker(children=sections,align='top',pad=0,sep=24),
            bbox_to_anchor=(.5,.025),bbox_transform=ax.transAxes,frameon=False)
        ax.add_artist(self.legend_box)
        return self.legend_box.get_window_extent(ax.figure.canvas.get_renderer()).height*72/ax.figure.dpi

    def _draw_graph(self, ax, positions, radii, scale):
        cfg = self.config['plot']
        renderer = ax.figure.canvas.get_renderer()
        label_boxes = []
        for u,v,data in self.network.edges(data=True):
            a,b = positions[u],positions[v]
            delta=b-a
            distance=np.linalg.norm(delta)
            if distance == 0:
                continue
            unit=delta/distance
            perp=np.array([-unit[1],unit[0]])
            color=cfg.get('edge_color','gray')
            start=a+unit*radii[u]; stop=b-unit*radii[v]
            curve=cfg.get('edge_style','line')=='curve'
            bend=.16*distance if curve else 0.
            t=np.linspace(0,1,32)
            line=start[None,:]+t[:,None]*(stop-start)+2*(t*(1-t))[:,None]*bend*perp
            ax.plot(*line.T,color=color,lw=cfg.get('edge_width',1),alpha=cfg.get('edge_alpha',.8),
                    ls='-',zorder=1)
            if cfg.get('show_edge_label',True):
                weight=data.get('weight',1)
                symbol=cfg.get('edge_label_symbol')
                label=symbol*int(weight) if symbol else f'{weight:g}'
                size=cfg.get('edge_label_size',8)
                tw,th,_=renderer.get_text_width_height_descent(label,FontProperties(size=size),False)
                half=np.array([tw,th])*36/ax.figure.dpi+2
                best=None
                for offset in (0,8,-8,16,-16):
                    for f in (.5,.38,.62,.25,.75):
                        xy=start+f*(stop-start)+(offset+2*f*(1-f)*bend)*perp
                        box=Bbox.from_extents(*(xy-half),*(xy+half))
                        penalty=sum(max(0,radii[n]+1-np.linalg.norm(np.maximum(
                            np.abs(xy-p)-half,0))) for n,p in positions.items())
                        penalty+=sum(10 for other in label_boxes if box.overlaps(other))
                        if best is None or penalty<best[0]:
                            best=(penalty,xy,box)
                        if penalty==0:
                            break
                    if best[0]==0: break
                _,xy,box=best
                label_boxes.append(box)
                ax.text(*xy,label,color=cfg.get('edge_label_color','black'),
                        size=size,ha='center',va='center',zorder=2,
                        bbox=dict(boxstyle='round,pad=.15',fc='white',ec='none'))
        for node,point in positions.items():
            self._draw_pie_node(ax,point,self.hap_composition.get(node,{}),
                                self.node_sizes[node],scale,cfg.get('node_alpha',.8))
            if cfg.get('show_node_label',False):
                ax.text(*point,str(node),fontsize=cfg.get('node_label_size',10),weight='semibold',
                        color=cfg.get('node_label_color','black'),ha='center',va='center',zorder=5)

    def plot(self):
        cfg=self.config['plot']
        width,height=72*cfg.get('width',10),72*cfg.get('height',8)
        fig=plt.figure(figsize=(width/72,height/72),facecolor='white')
        ax=fig.add_axes([0,0,1,1],xlim=(0,width),ylim=(0,height),aspect='equal')
        ax.set_axis_off()
        legend_height=self._legend(ax)
        bottom=height*.025+legend_height+(16 if legend_height else 6)
        title=cfg.get('title')
        top=height-(43 if title else 12)
        scale=self._area_scale()
        radii={n:np.sqrt(count*scale)/2 for n,count in self.node_sizes.items()}
        raw=self._get_layout_positions(cfg.get('layout','kamada_kawai'),cfg.get('seed'),cfg.get('k'))
        local=self._space_nodes(raw,radii,width,top-bottom)
        # Grow a crowded canvas instead of shrinking nodes.
        if local:
            xs=[p[0]+s*(radii[n]+6) for n,p in local.items() for s in (-1,1)]
            ys=[p[1]+s*(radii[n]+6) for n,p in local.items() for s in (-1,1)]
            padx=max(0,-min(xs),max(xs)-width)
            pady=max(0,-min(ys),max(ys)-(top-bottom))
            if padx or pady:
                width+=2*padx
                height+=2*pady; top+=2*pady
                local={n:p+[padx,pady] for n,p in local.items()}
                fig.set_size_inches(width/72,height/72)
                ax.set(xlim=(0,width),ylim=(0,height))
        if title: ax.text(12,height-8,title,va='top',fontsize=14,weight='semibold')
        self.node_positions=local
        pos={n:p+[0,bottom] for n,p in local.items()}
        self._draw_graph(ax,pos,radii,scale)
        if cfg.get('save_fig'):
            save_figure(fig,cfg['save_fig'],dpi=300,bbox_inches='tight',facecolor='white')
        plt.close(fig)
