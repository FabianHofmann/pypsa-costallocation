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
        snakemake = mock_snakemake('plot_payment_maps', nname='test-de10gf',
                                   method='ptpf', power='net')


n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions).set_index('name')

payments = xr.open_dataset(snakemake.input.payments)\
             .sum(['snapshot', 'source_carrier']).to_dataframe()

if not os.path.isdir(snakemake.output.folder):
    os.mkdir(snakemake.output.folder)


for col in payments:
    name = to_explanation[col]

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                           figsize=(5.5, 4))
    ax.spines['geo'].set_visible(False)
    regions.plot(column=payments[col], legend=True, ax=ax,
                 transform=ccrs.PlateCarree(), aspect='equal',
                 legend_kwds={'label': f'Payments for {name} [€]'})
    fig.canvas.draw()
    fig.tight_layout()
    fig.savefig(snakemake.output.folder + f'/{col}_total.png')


fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                        figsize=(5.5, 4))
ax.spines['geo'].set_visible(False)
nodal_payment = (n.buses_t.marginal_price * n.loads_t.p).sum()
regions.plot(column=nodal_payment, legend=True, ax=ax,
             transform=ccrs.PlateCarree(), aspect='equal',
             legend_kwds={'label': 'Payments [€]'})
fig.canvas.draw()
fig.tight_layout()
fig.savefig(snakemake.output.folder + '/electricity_total.png')
