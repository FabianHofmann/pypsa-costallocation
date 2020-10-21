#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 11:28:13 2020

@author: fabian
"""
import netallocation as ntl
import pypsa
import xarray as xr
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from numpy.random import rand
from pypsa.plot import compute_bbox_with_margins as compute_margins
import config
plt.rc('text', usetex=True)

if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('plot_graphical_abstract')


n = pypsa.Network(snakemake.input.network)
n.mremove('Link', n.links.index)
regions = gpd.read_file(snakemake.input.regions).set_index('name')

# %%

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8,5))
c = 'teal'
cmap = 'Greens'
sn = n.snapshots[451]
bus = 'DE0 10'
flowscale = 2e-2
busscale = 1e-5
bbox = dict(facecolor='w', alpha=1, boxstyle='round', pad=0.5, linewidth=0.5)

A = ntl.allocate_flow(n, sn)
P = A.peer_to_peer.sel(sink=bus).to_series()
F = A.peer_on_branch_to_peer.sum('source').sel(sink=bus).to_series()

lines_i = F[F.abs() > 0].Line.index
buses_i = n.lines.loc[lines_i, ['bus0', 'bus1']].stack().unique()


(x1, y1), (x2, y2) = compute_margins(0.2, n.buses.x[buses_i], n.buses.y[buses_i])
bounds = [x1, x2, y1, y2]

line_widths = pd.Series(1., lines_i).reindex(n.lines.index, fill_value=0)

gamma_line = (pd.Series(rand(lines_i.shape[0]), lines_i)
              .reindex(n.lines.index, fill_value=0))/2 + .5
gamma_bus = (pd.Series(rand(buses_i.shape[0]), buses_i)
             .reindex(n.buses.index, fill_value=0))/2 +.5

line_colors = pd.Series(c, lines_i).reindex(n.lines.index, fill_value='w')

n.plot(line_colors=gamma_line, line_cmap=cmap, line_widths=line_widths,
       bus_colors=gamma_bus, bus_cmap=cmap, bus_sizes=pd.Series(0.005, buses_i),
       geomap=False, ax=ax1, boundaries=bounds)
bbox.update({'edgecolor': 'darkgreen'})
ax1.set_title('Price intensity per\n network asset', bbox=bbox, color='darkgreen')
ax1.text(x2 , (y1+y2)/2, r'$+$', size=25, color='grey')


n.plot(flow=F * flowscale, bus_sizes=P * busscale, bus_colors=c,
       line_colors=line_colors, line_widths=line_widths,
       geomap=False, ax=ax2, boundaries=bounds,)
bbox.update({'edgecolor': c})
ax2.set_title('Dispatch \& Flow allocation\nto region', bbox=bbox, color=c)
ax2.text(x2 , (y1+y2)/2, r'$\rightarrow$', size=25, color='grey')
regions.loc[[bus]].plot(ax=ax2, aspect='equal', color='lightgrey')

n.plot(flow=F * flowscale, bus_sizes=P * busscale,
       line_colors=gamma_line, line_cmap=cmap, line_widths=line_widths,
       bus_colors=gamma_bus, bus_cmap=cmap,
       geomap=False, ax=ax3, boundaries=bounds)
bbox.update({'edgecolor': 'dimgrey'})
ax3.set_title('Cost allocation\n to region', bbox=bbox, color='k')
regions.loc[[bus]].plot(ax=ax3, aspect='equal', color='lightgrey')

fig.tight_layout(w_pad=0.5)
fig.savefig(snakemake.output[0], bbox_inches='tight', dpi=300)
