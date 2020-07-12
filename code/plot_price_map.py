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
        snakemake = mock_snakemake('plot_price_maps', clusters=50,
                                   method='ptpf', power='net')

n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions)
payments = xr.open_dataset(snakemake.input.payments)


demand = ntl.power_demand(n).rename(bus="payer")
bus_sizes = (payments).sum("snapshot").to_dataframe().stack().sort_index()
prices = (payments/demand).mean('snapshot').to_dataframe()
priceregions = regions.set_index('name').join(prices)

if not os.path.isdir(snakemake.output.folder):
    os.mkdir(snakemake.output.folder)


for col, name in to_explanation.items():

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                           figsize=(5, 4))
    ax.outline_patch.set_visible(False)
    priceregions.plot(column=col, legend=True, ax=ax,
                      transform=ccrs.PlateCarree(),
                      legend_kwds={'label': f'Average LMP for {name} [€/MWh]'})
    fig.tight_layout()
    fig.savefig(snakemake.output.folder + f'/{col}_average.png')
