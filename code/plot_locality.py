#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 17:47:01 2020

@author: fabian
"""
import netallocation as ntl
import pypsa
import xarray as xr
import geopandas as gpd
from config import to_explanation
import matplotlib.pyplot as plt
from helpers import combine_oneports, scfmt

if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('plot_locality', nname='test-de10bf',
                               method='ptpf', power='net')


n = pypsa.Network(snakemake.input.network)

coords = n.buses[['x', 'y']]
distance = (xr.DataArray(pypsa.geo.haversine(coords, coords),
                         coords=[n.buses.index, n.buses.index],
                         dims=['source', 'sink']))

cost = (combine_oneports(xr.open_dataset(snakemake.input.costs).sum('snapshot'))
        .drop('branch_investment_cost').drop_dims('branch')
        .assign(dist=distance))

df = (cost.stack(distance=['source', 'sink']).sortby('dist')
      .set_index(distance='dist').to_array().sum('variable').to_pandas().T)

cumshare = (df/df.sum()).sum(level=0).cumsum().mul(100) # in %
cumshare = cumshare[cumshare.iloc[20].sort_values(ascending=False).index]

# %%
fig, ax = plt.subplots(figsize=(5,5))
colors = n.carriers.color[cumshare.columns]
labels = n.carriers.nice_name[cumshare.columns]
cumshare.plot(ax=ax, color=colors, lw=2)
ax.set_xlabel('Distance [km]')
ax.set_ylabel('Cummulative Share of Payments [%]')
ax.legend(labels=labels, frameon=False, title=None)

fig.tight_layout()
fig.savefig(snakemake.output[0])
