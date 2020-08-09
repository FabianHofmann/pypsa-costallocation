#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 10:11:03 2020

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
        snakemake = mock_snakemake('plot_cost_maps', nname='test-de10gf',
                                   method='ptpf', power='net')

load_only = True

n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions)

payments = xr.open_dataset(snakemake.input.payments)

if load_only:
    payments = payments.sel(sink_carrier='Load', drop=True)
else:
    payments = payments.sum('sink_carrier')

payments = (payments).sum('snapshot').to_dataframe()
payments = payments.assign(lmp = (n.buses_t.marginal_price * n.loads_t.p).sum())
regions = regions.set_index('name').join(payments)

if not os.path.isdir(snakemake.output.folder):
    os.mkdir(snakemake.output.folder)


for col, name in to_explanation.items():

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                           figsize=(5.5, 4))
    ax.outline_patch.set_visible(False)
    regions.plot(column=col, legend=True, ax=ax,
                      transform=ccrs.PlateCarree(),
                      legend_kwds={'label': f'Total Payments for {name} [€]'})
    fig.canvas.draw()
    fig.tight_layout()
    fig.savefig(snakemake.output.folder + f'/{col}_total.png')


fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                        figsize=(5.5, 4))
ax.outline_patch.set_visible(False)
regions.plot(column='lmp', legend=True, ax=ax,
                  transform=ccrs.PlateCarree(),
                  legend_kwds={'label': f'Total Payments [€]'})
fig.canvas.draw()
fig.tight_layout()
fig.savefig(snakemake.output.folder + f'/electricity_total.png')
