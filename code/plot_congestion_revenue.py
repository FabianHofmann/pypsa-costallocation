#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 17:44:30 2020

@author: fabian
"""

import xarray as xr
import pypsa
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from helpers import scfmt


if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('plot_congestion_revenue', nname='de50')

ds = xr.load_dataset(snakemake.input.revenue).congestion_revenue
n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions)

b = ds.sum('payer')
line_widths = (b.sel(component='Line').to_series()
                .reindex(n.lines.index, fill_value=0)/b.sum().item()*10)
link_widths = (b.sel(component='Link').to_series()
                .reindex(n.links.index, fill_value=0)/b.sum().item()*10)

p = (ds.sum('branch').to_series().reindex(regions.index, fill_value=0))


fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()}, figsize=(5, 4))
ax.spines['geo'].set_visible(False)
n.plot(bus_sizes=0, line_widths=line_widths, link_widths=link_widths,
        ax=ax, boundaries=regions.total_bounds[[0,2,1,3]], geomap='10m')
regions.plot(column=p, transform=ccrs.PlateCarree(), aspect='equal', ax=ax,
              legend=True, legend_kwds={'label': ' [€]',
                                        'format': scfmt})
