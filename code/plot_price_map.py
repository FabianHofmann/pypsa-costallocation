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
from config import to_explanation
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_price_maps', nname='de10gf',
                                   method='ptpf', power='net')

n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions)

payments = xr.open_dataset(snakemake.input.payments)
demand = ntl.power_demand(n, per_carrier=True).rename(bus="sink")\
            .sel(carrier='Load', drop=True)

prices = (payments/demand).mean('snapshot').to_dataframe()
prices = prices.assign(lmp = n.buses_t.marginal_price.mean())
regions = regions.set_index('name').join(prices)
regions = regions.assign(demand = n.loads_t.p.mean())

if not os.path.isdir(snakemake.output.folder):
    os.mkdir(snakemake.output.folder)


for col, name in to_explanation.items():

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                           figsize=(5, 4))
    ax.outline_patch.set_visible(False)
    regions.plot(column=col, legend=True, ax=ax,
                      transform=ccrs.PlateCarree(),
                      legend_kwds={'label': f'Average LMP for \n {name} [€/MWh]'})
    fig.canvas.draw()
    fig.tight_layout()
    fig.savefig(snakemake.output.folder + f'/{col}_average.png')


fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                        figsize=(5, 4))
ax.outline_patch.set_visible(False)
regions.plot(column='lmp', legend=True, ax=ax,
                  transform=ccrs.PlateCarree(),
                  legend_kwds={'label': f'Average LMP [€/MWh]'})
fig.canvas.draw()
fig.tight_layout()
fig.savefig(snakemake.output.folder + f'/electricity_average.png')


fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                        figsize=(5, 4))
ax.outline_patch.set_visible(False)
regions.plot(column='demand', legend=True, ax=ax,
                  transform=ccrs.PlateCarree(),
                  legend_kwds={'label': f'Average Demand [MWh]'})
fig.canvas.draw()
fig.tight_layout()
fig.savefig(snakemake.output.folder + f'/demand_average.png')
