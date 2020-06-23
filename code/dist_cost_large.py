#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:13:56 2020

@author: fabian
"""

import pypsa
import netallocation as ntl
from  netallocation.cost import (nodal_demand_cost, nodal_production_revenue)
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from config import to_symbol, color
plt.rc('text', usetex=True)


n = pypsa.Network('../elec_s_50_ec_lvopt_Co2L-3H.nc')
n.calculate_dependent_values()
n.snapshot_weightings[:] = 1
n.generators['p_nom_extendable'] = True
n.lines['s_nom_min'] = 0
# n.lines['s_nom_max'] = 6000
# n.generators['p_nom_max'] = np.inf
n.mremove('StorageUnit', n.storage_units.index)
n.mremove('Link', n.links.index)
n.mremove('Line', n.lines.query('x_pu == inf').index)
n.remove('Bus', 'DE0 0')

T = 10
n.set_snapshots(n.snapshots[:T])
n.global_constraints.loc['CO2Limit', 'constant'] *= T/8700
n.lopf(pyomo=False, keep_shadowprices=True, solver_name='gurobi')
tag = '_large'


# %%
ds = ntl.allocate_flow(n, method='ptpf', aggregated=False, #q=0
                       )
pr = nodal_production_revenue(n).rename(bus='payer')
dc = nodal_demand_cost(n).rename(bus='payer')
ca = ntl.allocate_cost(n, method=ds)
payments = ca.sum([d for d in ca.dims if d not in ['payer', 'snapshot']])
dc = nodal_demand_cost(n).rename(bus='payer')


# %% Network plot

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,5))
bus_sizes = n.generators.groupby(['bus', 'carrier']).p_nom_opt.sum()
plot = n.plot(bus_sizes=bus_sizes/1e5, line_widths=n.lines.s_nom_opt/1e3,
        link_widths=n.links.p_nom_opt/1e3, ax=ax)
ax.legend(*ntl.plot.handles_labels_for(n.carriers.color), loc='center left',
          bbox_to_anchor=(1, 0.5), ncol=1)
fig.canvas.draw(); fig.tight_layout()
fig.savefig(f'../figures/network{tag}.png')


# %% Allocation plot


fig, axes = plt.subplots(3,1, figsize=(12,7), sharex=True)

for sn, ax in zip(n.snapshots, axes.ravel()):
    ca_df = payments.sel(snapshot=sn, drop=True).to_dataframe()
    dc_df = dc.sel(snapshot=sn, drop=True).to_dataframe()
    d = dict(stacked=True, zorder=3, width=.3, legend=False, ax=ax)

    ca_df.plot.bar(position=1.1, color=color[ca_df.columns], **d)
    dc_df.plot.bar(position=0, color=color[dc_df.columns], **d)

    ax.set_xlim(left=-.5)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(5,2))
    ax.set_xlabel('')
    ax.set_title(f'$t = %i$'%(sn.hour+1))
    ax.grid(axis='y', zorder=1, linestyle='dashed')

handles, labels = ax.get_legend_handles_labels()
labels = [to_symbol[l] for l in labels]
fig.legend(handles[:-1], labels[:-1], ncol=4, frameon=False, loc='lower center',
            bbox_to_anchor=(0.3, 1), fontsize='large', title='Flow Based (left bars)')
fig.legend(handles[-1:], labels[-1:], ncol=1, frameon=False, loc='lower center',
            bbox_to_anchor=(0.7, 1), fontsize='large', title='LMP Based (right bars)')
fig.tight_layout()
fig.savefig(f'../figures/compare_allocation{tag}.png', bbox_inches='tight')

#%% LMP

demand = ntl.power_demand(n).rename(bus='payer')
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(5,4))
bus_sizes = (payments).sum('snapshot').to_dataframe().stack().sort_index()
n.plot(ax=ax, bus_sizes=bus_sizes/bus_sizes.sum()*3, bus_colors=color)
handles, labels = ntl.plot.handles_labels_for(color[:-1])
labels = ['$\sum_t' + to_symbol[l][1:] for l in labels]
ax.legend(handles, labels, loc='upper left')
fig.canvas.draw(); fig.tight_layout()
fig.savefig(f'../figures/nodal_payments{tag}.png')


