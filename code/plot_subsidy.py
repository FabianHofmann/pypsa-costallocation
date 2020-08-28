#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 16:59:13 2020

@author: fabian
"""
import pypsa
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import cartopy.crs as ccrs
import netallocation as ntl
from matplotlib.pyplot import Line2D
from pypsa.plot import projected_area_factor
from netallocation.plot_helpers import (make_handler_map_to_scale_circles_as_in,
                                        make_legend_circles_for)

if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('plot_subsidy', nname='test-de10bf')


n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions).set_index('name')


oneports = pd.concat([n.df(c).groupby(['bus', 'carrier']).mu_lower_p_nom.sum()
                      for c in ['Generator', 'StorageUnit']]).sort_index()

line_widths  = n.lines.mu_lower_s_nom
link_widths = n.links.mu_lower_p_nom
branch_sum = line_widths.sum() + link_widths.sum()

bus_scale = 2
branch_scale = 100

fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                        figsize=(5, 5))
ax.spines['geo'].set_visible(False)

n.plot(bus_sizes=oneports/oneports.sum() * bus_scale,
       line_widths = line_widths / branch_sum * branch_scale,
       link_widths = link_widths / branch_sum * branch_scale,
       ax=ax,
       geomap='10m',
       boundaries=regions.total_bounds[[0,2,1,3]],
       )
regions.plot(ax=ax, transform=ccrs.PlateCarree(), aspect='equal',
             color='white', lw=0.3, edgecolor='grey')


fig.legend(
    *ntl.plot.handles_labels_for(n.carriers.set_index('nice_name').color),
    loc="upper left",
    bbox_to_anchor=(1, 1),
    title='Carrier',
    frameon=False,
)


# legend generator capacities
reference_caps = [50e3, 10e3]
scale = oneports.sum() / bus_scale / projected_area_factor(ax)**2
handles = make_legend_circles_for(reference_caps, scale=scale,
                                  facecolor="w", edgecolor='grey',
                                  alpha=.5)
labels = ["%ik €"%(s / 1e3) for s in reference_caps]
handler_map = make_handler_map_to_scale_circles_as_in(ax)
legend = fig.legend(handles, labels,
                loc="upper left", bbox_to_anchor=(1., 0.4),
                frameon=False,
                handler_map=handler_map)
fig.add_artist(legend)

# legend transmission capacitues
handles, labels = [], []
reference_caps = [10e3, 5e3]
handles = [Line2D([0], [0], color='grey', linewidth=s * branch_scale / branch_sum)
           for s in reference_caps]
labels = ['%ik €'%(s / 1e3) for s in reference_caps]

legend = fig.legend(handles, labels,
                    loc="lower left", bbox_to_anchor=(1, .1),
                    frameon=False,
                    )
fig.artists.append(legend)


fig.canvas.draw()
fig.tight_layout()
fig.savefig(snakemake.output[0], bbox_inches='tight')
