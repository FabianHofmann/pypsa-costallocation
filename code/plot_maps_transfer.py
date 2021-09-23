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
from helpers import combine_oneports, scfmt, load

if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('plot_price_maps', nname='de50')

plt.rc('figure', dpi=300)

n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions).set_index('name')
cost = (combine_oneports(xr.open_dataset(snakemake.input.costs).sum('snapshot'))
        .set_index({'branch': ['component', 'branch_i']})
        .transpose('sink', 'source', 'source_carrier', 'branch'))


if 'price' in snakemake.output[0]:
    if any(n.snapshot_weightings.objective != 1):
        w = ntl.cost.snapshot_weightings(n)
        cost['one_port_operational_cost'] /= w
        cost['co2_cost'] /= w
    payment_type = 'Average Price for \n'
    unit = '€/MWh'
    demand = load(n).sum('snapshot')
    cost = cost / demand
    fmt = None

else:
    payment_type = 'Allocated'
    unit = '€'
    fmt = scfmt

payer = cost.sum(['source', 'branch']).reindex(sink=regions.index).fillna(0)
receiver = cost.sum('sink').reindex(source=regions.index).fillna(0)

res = '50m' if 'test' in snakemake.output[0] else '10m'
nplot_kwargs = dict(bus_colors='rosybrown', geomap=res,
                    line_widths=0, link_widths=0, bus_alpha=0.8,
                    boundaries=regions.total_bounds[[0,2,1,3]])
rplot_kwargs = dict(legend=True, transform=ccrs.PlateCarree(), aspect='equal')

os.makedirs(snakemake.output.folder + '/by_carrier', exist_ok=True)

for var in cost:
    p = payer[var].to_pandas()
    r = receiver[var].to_pandas()
    varname = to_explanation[var]
    for carrier in cost.source_carrier.data:
        if 'branch' in var: continue
        if (p[carrier].sum() <= 0.001): continue

        expand = varname.replace('Production & Storage ', '')
        ncarrier = n.carriers.nice_name[carrier]

        fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                                figsize=(5, 4))
        ax.spines['geo'].set_visible(False)

        n.plot(bus_sizes=r[carrier]/r[carrier].sum()*0.7, ax=ax, **nplot_kwargs)
        regions.plot(column=p[carrier], ax=ax, **rplot_kwargs,
                     legend_kwds={'label': f'{payment_type} {ncarrier} {expand} [{unit}]',
                                  'format': fmt})

        fig.canvas.draw(), fig.tight_layout()
        fig.savefig(snakemake.output.folder + f'/by_carrier/{carrier}_{var}.png')


    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                           figsize=(5, 4))
    ax.spines['geo'].set_visible(False)

    if 'branch' in var:
        normed = r/r.sum()
        n.plot(bus_sizes=0, ax=ax,
               line_widths=normed['Line'] * 80,
               link_widths=normed['Link'] * 80,
               boundaries=regions.total_bounds[[0,2,1,3]], geomap=res)
    else:
        n.plot(bus_sizes=r.sum(1)/p.sum().sum(), ax=ax, **nplot_kwargs)
    regions.plot(column=p if 'branch' in var else p.sum(1),
                 ax=ax, **rplot_kwargs,
                 legend_kwds={'label': f'{payment_type} {varname} [{unit}]',
                              'format': fmt})

    fig.canvas.draw(), fig.tight_layout()
    fig.savefig(snakemake.output.folder + f'/{var}.png')

