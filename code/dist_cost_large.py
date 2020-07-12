#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:13:56 2020

@author: fabian
"""

import pypsa
import netallocation as ntl
from netallocation.cost import nodal_demand_cost, nodal_production_revenue
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from config import to_symbol, color, to_total_symbol, to_explanation
import geopandas as gpd

plt.rc("text", usetex=True)
plt.rc('font', family='serif')


n = pypsa.Network('../solved_germany_50.nc')
n.set_snapshots(n.snapshots[:10])
tag = "_large"


# %% total cost
plt.rc("text", usetex=False)
ca.sum().to_array().to_series().reindex_like(color).dropna()\
    .rename(to_total_symbol).rename('')\
    .plot.pie(colors=color, explode=[.01]*len(ca), autopct='%1.0f%%')
plt.rc("text", usetex=True)




#%% LMP
demand = ntl.power_demand(n).rename(bus="payer")
bus_sizes = (payments).sum("snapshot").to_dataframe().stack().sort_index()
prices = (payments/demand).mean('snapshot').to_dataframe()
regions = gpd.read_file('../regions_onshore_elec_s_50.geojson')
priceregions = regions.set_index('name').join(prices)


for col, name in to_explanation.items():

    fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                           figsize=(5, 4))
    ax.outline_patch.set_visible(False)
    priceregions.plot(column=col, legend=True, ax=ax,
                      transform=ccrs.PlateCarree(),
                      legend_kwds={'label': f'Average LMP for {name} [€/MWh]'})
    fig.tight_layout()
    fig.savefig(f'../figures/{col}_average.png')

# %%

fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                       figsize=(5, 4))
n.plot(ax=ax, bus_sizes=bus_sizes / bus_sizes.sum() * 3, bus_colors=color)
handles, labels = ntl.plot.handles_labels_for(color[bus_sizes.index.unique(1)])
labels = ["$\sum_t" + to_symbol[l][1:] for l in labels]
ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1, 1), ncol=1, frameon=False)
fig.canvas.draw()
fig.tight_layout()
fig.savefig(f"../figures/nodal_payments{tag}.png")
