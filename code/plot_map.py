#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 23:09:45 2020

@author: fabian
"""

import netallocation as ntl
import pypsa
import xarray as xr
import geopandas as gpd
import cartopy.crs as ccrs
from config import to_explanation, source_carrier
import matplotlib.pyplot as plt
import os

if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_average_demand', nname='de50')

n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions).set_index('name')

if 'price' in snakemake.output[0]:
    data = n.buses_t.marginal_price.mean()
    label = 'Average LMP [€/MWh]'
elif 'demand' in snakemake.output[0]:
    data = n.loads_t.p.mean()
    label = 'Average Demand [MW]'

data = data.reindex(regions.index).fillna(0)

fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                        figsize=(5, 4))
ax.spines['geo'].set_visible(False)
regions.plot(column=data, legend=True, ax=ax,
             transform=ccrs.PlateCarree(), aspect='equal',
             legend_kwds={'label': label})
fig.canvas.draw()
fig.tight_layout()
fig.savefig(snakemake.output[0])

