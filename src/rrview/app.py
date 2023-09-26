import argparse
import xarray as xr
import numpy as np
import holoviews as hv
import hvplot.xarray
import panel as pn
import param
import geoviews as gv
from pathlib import Path
import geonetworkx as gnx
import networkx as nx
from holoviews import opts
import warnings

from bokeh.models.formatters import NumeralTickFormatter
from bokeh.models import HoverTool

warnings.filterwarnings('ignore')
hv.extension('bokeh')

template = pn.template.BootstrapTemplate(title='River-regulation viewer')
get_downstream = lambda node: list(G.successors(node))
get_upstreams = lambda node: list(G.predecessors(node))

parser = argparse.ArgumentParser(description='River regulation viewer')
parser.add_argument('rrds_fn', type=str, help='Path to data ')
parser.add_argument('-nlnk', '--network-links', type=str, help='Path to reservoir network')
parser.add_argument('-npts', '--network-points', type=str, help='Path to reservoir network points')

args = parser.parse_args()

def setup_data(rrds_fn, nlnk_fn, npts_fn):
    global ds
    global G
    global nodes
    global default_current_node
    global default_downstream_node
    global default_upstream_nodes

    ds = xr.open_dataset(rrds_fn)
    G = gnx.read_geofiles(npts_fn, nlnk_fn, directed=True)
    nodes = ds['node'].values.tolist()
    default_current_node = nodes[0]
    default_downstream_node = get_downstream(default_current_node)
    default_upstream_nodes = get_upstreams(default_current_node)

setup_data(args.rrds_fn, args.network_links, args.network_points)

def plot_node(node):
    default_opts = dict(
        height=400, width=800, 
        yformatter=NumeralTickFormatter(format='0.00a'), 
        ylabel='Flow (m3/d)', xlabel='Date',
    )
    default_curve_opts = dict(
        interpolation='steps-pre'
    )
    
    hover_cols = ['inflow', 'outflow', 'theoretical_natural_runoff', 'storage_change', 'obs_inflow']
    tooltips = [
        ('Inflow', '@inflow{0.00a}'),
        ('Observed Inflow', '@obs_inflow{0.00a}'),
        # ('Outflow', '@outflow{0.00a}'),
        ('TNR', '@theoretical_natural_runoff{0.00a}'),
        ('Storage Change', '@storage_change{0.00a}')
    ]
    hover = HoverTool(tooltips=tooltips, mode='vline')

    inflow_plot = ds.sel(node=node) \
        .hvplot(kind="line", x='time', y='inflow', label='Inflow', hover_cols=hover_cols) \
        .opts(**default_opts) \
        .opts(**default_curve_opts) \
        .opts(tools=[hover])
    
    outflow_plot = ds.sel(node=node) \
        .hvplot(kind="line", x='time', y='outflow', label='Outflow') \
        .opts(**default_opts) \
        .opts(**default_curve_opts) \
        .opts(muted=True)
    
    tnr_plot = ds.sel(node=node) \
        .hvplot(kind="line", x='time', y='theoretical_natural_runoff', label='Theoretical Natural Runoff (TNR)') \
        .opts(**default_opts) \
        .opts(**default_curve_opts) \
        .opts(color='green', line_dash='2 2')
    
    obs_inflow_plot = ds.sel(node=node) \
        .hvplot(kind="line", x="time", y='obs_inflow', label='Observed Inflow') \
        .opts(**default_opts) \
        .opts(**default_curve_opts) \
        .opts(color='k', line_dash='2 2', alpha=0.8)
    
    ds['zeros'] = 0
    dels_plot = ds.sel(node=node) \
        .hvplot(kind="area", x='time', y='storage_change', y2='zeros', label='Storage Change') \
        .opts(**default_opts) \
        .opts(fill_alpha=0.2, color='gray', muted_fill_alpha=0.05, muted=True, muted_line_alpha=0.3)
    

    return inflow_plot * outflow_plot * tnr_plot * obs_inflow_plot * dels_plot


class RegulationViewer(param.Parameterized):
    current_node = param.Selector(objects=nodes, precedence=1)
    downstream_node = param.Number(default=default_downstream_node[0] if len(default_downstream_node) else None, allow_None=True, precedence=2)
    upstream_nodes = param.Selector(default=default_upstream_nodes[0] if len(default_upstream_nodes) else None, objects=default_upstream_nodes, allow_None=True, precedence=0)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.current_hv = hv.Curve([])
        self.upstream_hv = hv.Curve([])
        self.downstream_hv = hv.Curve([])

        self._update_upstream_nodes()
        self._update_downstream_node()

        self._plot_current_node()
        self._plot_upstream_node()
        self._plot_downstream_node()

    @param.depends('current_node', watch=True)
    def _update_downstream_node(self):
        downstreams = get_downstream(self.current_node)
        self.downstream_node = downstreams[0] if len(downstreams)>0 else None
    
    @param.depends('current_node', watch=True)
    def _update_upstream_nodes(self):
        self.param.upstream_nodes.objects = get_upstreams(self.current_node)
        if len(self.param.upstream_nodes.objects)>0:
            self.upstream_nodes = self.param.upstream_nodes.objects[0]
        else:
            self.upstream_nodes = None

    @param.depends('current_node', watch=True)
    def _plot_current_node(self):
        self.current_hv = plot_node(self.current_node)

    @param.depends('upstream_nodes', watch=True)
    def _plot_upstream_node(self):
        if self.upstream_nodes is not None:
            self.upstream_hv = plot_node(self.upstream_nodes)
        else:
            self.upstream_hv = hv.Curve([])

    @param.depends('downstream_node', watch=True)
    def _plot_downstream_node(self):
        if self.downstream_node is not None:
            self.downstream_hv = plot_node(self.downstream_node)
        else:
            self.downstream_hv = hv.Curve([])

    def _style_current_node(self):
        return self.current_hv.opts(title=f'Current Node: {self.current_node} - {G.nodes[self.current_node]["name"].replace("_", " ")}')
    
    def _style_upstream_node(self):
        if self.upstream_nodes:
            title = f'Upstream Node: {self.upstream_nodes} - {G.nodes[self.upstream_nodes]["name"].replace("_", " ")}'
        else:
            title = f'No upstream nodes'
        
        return self.upstream_hv.opts(title=title)
    
    def _style_downstream_node(self):
        if self.downstream_node:
            title = f'Downstream Node: {self.downstream_node} - {G.nodes[self.downstream_node]["name"].replace("_", " ")}'
        else:
            title = f'No downstream nodes'
        return self.downstream_hv.opts(title=title)

    @param.depends('current_node')
    def _info_plot(self):
        upstream_nodes_str = ", ".join(map(str, self.param.upstream_nodes.objects))
        downstream_nodes_str = str(self.downstream_node)
        
        return pn.pane.Markdown(f"""
                    *Upstream node(s)*: {upstream_nodes_str}
                    **Current node**: {self.current_node}
                    *Downstream node(s)*: {downstream_nodes_str}
                    """
                )
    
    @property
    def params_pane(self):
        return pn.Param(self.param, widgets={'downstream_node': {'type': pn.widgets.IntInput, 'disabled': True}})


    @param.depends('current_node', 'upstream_nodes', 'current_hv', 'upstream_hv')
    def view(self):
        return pn.Row(
            pn.Column(
                self._info_plot, self.params_pane
            ),
            pn.Column(
                self._style_upstream_node,
                self._style_current_node,
                self._style_downstream_node
            )
        )

def setup_dashboard(rv):
    global template

    template.sidebar.append(
        pn.Column(rv._info_plot, rv.params_pane)
    )
    template.main.append(
        pn.Column(
            rv._style_upstream_node,
            rv._style_current_node,
            rv._style_downstream_node
        )
    )
    return template.servable();


rv = RegulationViewer()
setup_dashboard(rv)