#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 17:44:30 2020

@author: fabian
"""



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
