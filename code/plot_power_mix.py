#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 10:48:13 2020

@author: fabian
"""
import pypsa
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import netallocation as ntl
import config
import xarray as xr

if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('plot_power_mix', nname='de50')

plt.rc('figure', dpi=300)

n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions)
power_mix = xr.open_dataarray(snakemake.input.power_mix)

bus_sizes = (power_mix / power_mix.sum('source_carrier')).to_series() /5


fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                         figsize=(5, 5))
# optimal
plot = n.plot(bus_sizes=bus_sizes.fillna(0)*.7, line_widths=0, link_widths=0,
              ax=ax, geomap='10m',
              boundaries=regions.total_bounds[[0,2,1,3]] )

regions.plot(ax=ax, transform=ccrs.PlateCarree(), aspect='equal',
             color='white', edgecolor='lightgrey', lw=0.5)

fig.legend(
    *ntl.plot.handles_labels_for(n.carriers.set_index('nice_name').color),
    loc="lower center",
    bbox_to_anchor=(.5, 1),
    title='Carrier',
    frameon=False,
    ncol=2
)



fig.canvas.draw()
fig.tight_layout()
fig.savefig(snakemake.output[0], bbox_inches='tight')