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
import netallocation as ntl

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_price_maps', nname='test-de10gf',
                                   method='ptpf', power='net')

if not os.path.isdir(snakemake.output.folder):
    os.mkdir(snakemake.output.folder)

n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions).set_index('name')

payments = xr.open_dataset(snakemake.input.payments).sum('source_carrier')

w = ntl.cost.snapshot_weightings(n)
payments['one_port_operational_cost'] /= w
payments['co2_cost'] /= w

demand = ntl.power_demand(n, per_carrier=True).sel(carrier='Load', drop=True)\
            .rename(bus="sink")
prices = (payments/demand).mean('snapshot').to_dataframe()


for col in prices:
    name = to_explanation[col]

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                           figsize=(5, 4))
    ax.spines['geo'].set_visible(False)
    regions.plot(column=prices[col], legend=True, ax=ax,
                 transform=ccrs.PlateCarree(), aspect='equal',
                 legend_kwds={'label': f'Average LMP for \n {name} [€/MWh]'})
    fig.canvas.draw()
    fig.tight_layout()
    fig.savefig(snakemake.output.folder + f'/{col}_average.png')


fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                        figsize=(5, 4))
ax.spines['geo'].set_visible(False)
lmp = n.buses_t.marginal_price.mean()
regions.plot(column=lmp, legend=True, ax=ax,
             transform=ccrs.PlateCarree(), aspect='equal',
             legend_kwds={'label': 'Average LMP [€/MWh]'})
fig.canvas.draw()
fig.tight_layout()
fig.savefig(snakemake.output.folder + '/electricity_average.png')


fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                        figsize=(5, 4))
ax.spines['geo'].set_visible(False)
demand = n.loads_t.p.mean()
regions.plot(column=demand, legend=True, ax=ax,
             transform=ccrs.PlateCarree(), aspect='equal',
             legend_kwds={'label': 'Average Demand [MWh]'})
fig.canvas.draw()
fig.tight_layout()
fig.savefig(snakemake.output.folder + '/demand_average.png')
