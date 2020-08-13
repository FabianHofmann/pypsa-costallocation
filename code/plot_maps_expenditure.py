#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 10:36:46 2020

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
        snakemake = mock_snakemake('plot_expenditure_maps', nname='test-de10gf',
                                   method='ptpf', power='net')

n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions)

expenditures = xr.open_dataset(snakemake.input.costs).sum(['sink', 'snapshot'])\
                 .sel(sink_carrier='Load', drop=True) \
                 .rename(source='bus', source_carrier='carrier')\
                 .set_index(branch=['component', 'branch_i'])

by_bus_carrier = (expenditures.drop('branch_investment_cost')
                  .drop_dims('branch').to_dataframe().unstack('carrier'))
by_bus = by_bus_carrier.sum(level=0, axis=1)
by_branch = expenditures.branch_investment_cost.to_series()

# regions = regions.assign(production = production.mean())

if not os.path.isdir(snakemake.output.folder):
    os.mkdir(snakemake.output.folder)

for col in by_bus:
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                            figsize=(5, 4))
    ax.spines['geo'].set_visible(False)
    regions.plot(column=by_bus[col], legend=True, ax=ax,
                 transform=ccrs.PlateCarree(), aspect='equal',
                 legend_kwds={'label': f'Total {to_explanation[col]} [€]'})
    fig.canvas.draw()
    fig.tight_layout()
    fig.savefig(snakemake.output.folder + f'/{col}_total.png')


if not os.path.isdir(snakemake.output.folder + '/by_carrier'):
    os.mkdir(snakemake.output.folder + '/by_carrier')

for col in by_bus_carrier:
    carrier = n.carriers.nice_names[col[1]]
    expend = to_explanation[col[0]]
    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                            figsize=(5, 4))
    ax.spines['geo'].set_visible(False)
    regions.plot(column=by_bus_carrier[col], legend=True, ax=ax,
                 transform=ccrs.PlateCarree(), aspect='equal',
                 legend_kwds={'label': f'{carrier} {expend} [€]'})
    fig.canvas.draw()
    fig.tight_layout()
    fig.savefig(snakemake.output.folder + f'/by_carrier/{col[1]}_{col[0]}.png')



# fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
#                         figsize=(5, 4))
# ax.outline_patch.set_visible(False)
# regions.plot(column='lmp', legend=True, ax=ax,
#                   transform=ccrs.PlateCarree(),
#                   legend_kwds={'label': f'Average LMP [€/MWh]'})
# fig.canvas.draw()
# fig.tight_layout()
# fig.savefig(snakemake.output.folder + f'/electricity_average.png')
