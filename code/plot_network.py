#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 10:39:11 2020

@author: fabian
"""

import pypsa
import config
import matplotlib.pyplot as plt
import netallocation as ntl
import cartopy.crs as ccrs
import geopandas as gpd
from matplotlib.pyplot import Line2D
from netallocation.plot_helpers import (make_handler_map_to_scale_circles_as_in,
                                        make_legend_circles_for)
from pypsa.plot import projected_area_factor

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_network', nname='de50')


plt.rc('axes', titlesize='medium')

n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions)

bus_scale = 1.5e-5
branch_scale = 5e-4

fig, axes = plt.subplots(1, 2, subplot_kw={"projection": ccrs.EqualEarth()},
                         figsize=(8, 4))
# existent
bus_sizes = n.generators.groupby(["bus", "carrier"]).p_nom_min.sum()
plot = n.plot(
    bus_sizes=bus_sizes * bus_scale,
    line_widths=n.lines.s_nom_min * branch_scale,
    link_widths=n.links.p_nom_min * branch_scale,
    ax=axes[0],
    geomap='10m',
    title='Lower Capacity Bounds',
    boundaries=regions.total_bounds[[0,2,1,3]]
)
regions.plot(ax=axes[0], transform=ccrs.PlateCarree(), aspect='equal', facecolor='white',
             edgecolor='blue', linewidth=0.1)

# expanded
bus_sizes = n.generators.eval('p_nom_ext = p_nom_opt - p_nom_min')\
             .groupby(["bus", "carrier"]).p_nom_ext.sum()
plot = n.plot(
    bus_sizes=bus_sizes * bus_scale,
    line_widths=n.lines.eval('s_nom_opt - s_nom_min') * branch_scale,
    link_widths=n.links.eval('p_nom_opt - p_nom_min') * branch_scale,
    ax=axes[1],
    geomap='10m',
    title='Capacity Expansion',
    boundaries=regions.total_bounds[[0,2,1,3]]
)
regions.plot(ax=axes[1], transform=ccrs.PlateCarree(), aspect='equal', facecolor='white',
             edgecolor='blue', linewidth=0.1)

fig.legend(
    *ntl.plot.handles_labels_for(n.carriers.set_index('nice_name').color),
    loc="upper left",
    bbox_to_anchor=(0,0),
    ncol=2,
    title='Generation Type',
    frameon=False,
)

# legend generator capacities
reference_caps = [10e3, 5e3, 1e3]
scale = 1 / bus_scale / projected_area_factor(axes[0])**2
handles = make_legend_circles_for(reference_caps, scale=scale,
                                  facecolor="w", edgecolor='grey',
                                  alpha=.5)
labels = ["%i GW"%(s / 1e3) for s in reference_caps]
handler_map = make_handler_map_to_scale_circles_as_in(axes[0])
legend = fig.legend(handles, labels,
                loc="upper left", bbox_to_anchor=(.5, 0),
                frameon=False,  # edgecolor='w',
                title='Generation Capacity',
                handler_map=handler_map)
fig.add_artist(legend)



# legend AC / DC
handles = [Line2D([0], [0], color=c, linewidth=5) for c in
           ['rosybrown', 'darkseagreen']]
labels = ['AC', 'DC']

legend = fig.legend(handles, labels,
                    loc="upper left", bbox_to_anchor=(0.8, 0),
                    frameon=False,
                    title='Transmission Type')
fig.artists.append(legend)


# legend transmission capacitues
handles, labels = [], []
reference_caps = [10,5]
handles = [Line2D([0], [0], color='grey', alpha=0.5, linewidth=s*1e3*branch_scale)
           for s in reference_caps]
labels = ['%i GW'%s for s in reference_caps]

legend = fig.legend(handles, labels,
                    loc="upper left", bbox_to_anchor=(.8, -.15),
                    frameon=False,
                    title='Transmission Capacity')
fig.artists.append(legend)



fig.canvas.draw()
fig.tight_layout()
fig.savefig(snakemake.output[0], bbox_inches='tight')



