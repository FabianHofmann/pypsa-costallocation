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
from config import to_explanation, source_dims_r
import matplotlib.pyplot as plt
import os
from helpers import combine_oneports, fmt

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_expenditure_maps', nname='test-de10bf',
                                   method='ptpf', power='net')

n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions).set_index('name')

cost = combine_oneports(xr.open_dataset(snakemake.input.costs).sum('snapshot'))
cost = cost.set_index({'branch': ['component', 'branch_i']})
cost = cost.transpose('sink', 'source', 'source_carrier', 'branch')

payer = cost.sum('source')
receiver = cost.sum('sink')


if not os.path.isdir(snakemake.output.folder):
    os.mkdir(snakemake.output.folder)


for var in cost:
    if var == 'branch_investment_cost': continue
    p = payer[var].to_pandas().reindex(regions.index).fillna(0)
    r = receiver[var].to_pandas().reindex(regions.index).fillna(0)
    varname = to_explanation[var]
    for carrier in cost.source_carrier.data:
        if p[carrier].sum() <= 0.001: continue

        expend = varname.replace('Production & Storage ', '')
        ncarrier = n.carriers.nice_name[carrier]

        fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                                figsize=(5, 4))
        ax.spines['geo'].set_visible(False)

        n.plot(bus_sizes=r[carrier]/r[carrier].sum(),
               bus_colors='rosybrown',
               line_widths=0, link_widths=0, ax=ax,
               boundaries=regions.total_bounds[[0,2,1,3]])
        regions.plot(column=p[carrier],
                     legend=True, ax=ax,
                     transform=ccrs.PlateCarree(), aspect='equal',
                     legend_kwds={'label': f'Allocated {ncarrier} {expend} [€]',
                                  'format': fmt})

        fig.canvas.draw()
        fig.tight_layout()
        fig.savefig(snakemake.output.folder + f'/by_carrier/{carrier}_{var}.png')


    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                            figsize=(5, 4))
    ax.spines['geo'].set_visible(False)

    n.plot(bus_sizes=p.sum(1)/p.sum().sum(),
           bus_colors='rosybrown',
           line_widths=0, link_widths=0, ax=ax,
           boundaries=regions.total_bounds[[0,2,1,3]])
    regions.plot(column=r.sum(1),
                 legend=True, ax=ax,
                 transform=ccrs.PlateCarree(), aspect='equal',
                 legend_kwds={'label': f'Allocated {varname} [€]', 'format': fmt})

    fig.canvas.draw()
    fig.tight_layout()
    fig.savefig(snakemake.output.folder + f'/{var}_total.png')

# %%
expanded_i = n.branches().query('s_nom_opt - s_nom_min > 1 | '
                                'p_nom_opt - p_nom_min > 1').index

b = cost.branch_investment_cost.sel(branch=expanded_i).sum('sink')
line_widths = (b.sel(component='Line').to_series()
               .reindex(n.lines.index, fill_value=0)/b.sum().item()*10)
link_widths = (b.sel(component='Link').to_series()
               .reindex(n.links.index, fill_value=0)/b.sum().item()*10)

p = (cost.branch_investment_cost.sel(branch=expanded_i).sum('branch')
    .to_series().reindex(regions.index, fill_value=0))


fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()}, figsize=(5, 4))
ax.spines['geo'].set_visible(False)
n.plot(bus_sizes=0, line_widths=line_widths, link_widths=link_widths,
       ax=ax, boundaries=regions.total_bounds[[0,2,1,3]], geomap='10m')
regions.plot(column=p, transform=ccrs.PlateCarree(), aspect='equal', ax=ax,
             legend=True, legend_kwds={'label': 'Payment to Expanded Transmission [€]',
                                       'format': fmt})
