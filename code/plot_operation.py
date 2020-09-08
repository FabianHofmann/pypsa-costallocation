#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 15:43:11 2020

@author: fabian
"""
import netallocation as ntl
import pypsa
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import geopandas as gpd
import config

if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('plot_operation_high_expenditure', nname='test-de10bf')


n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions).set_index('name')

sn = (n.loads_t.p * n.buses_t.marginal_price).sum(1).idxmax()

# %%
fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                        figsize=(5, 4))

prod = ntl.power_production(n, sn, per_carrier=True).to_pandas().stack()

n.plot(flow=sn, line_widths=0.01, link_widths=.01, bus_sizes=prod/1e5, ax=ax,
       boundaries=regions.total_bounds[[0,2,1,3]], geomap='10m')
regions.plot(column=n.loads_t.p.loc[sn], ax=ax, legend=True,
             transform=ccrs.PlateCarree(), aspect='equal', cmap='Purples',
             alpha=0.5, legend_kwds={'alpha':.5, 'label': 'Load [MW]'}, vmin=0)

fig.canvas.draw()
fig.tight_layout()
fig.savefig(snakemake.output[0], bbox_inches='tight')
