#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:13:56 2020

@author: fabian
"""

import pypsa
import pandas as pd
import numpy as np
import netallocation as ntl
from  netallocation.cost import (nodal_co2_cost, nodal_demand_cost,
                                 nodal_production_revenue, congestion_revenue)
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
plt.rc('text', usetex=True)


to_symbol_dict = dict(one_port_investment_cost = '$\mathcal{C}^{G}_{n,t}$',
                      branch_investment_cost = '$\mathcal{C}^{F}_{n,t}$',
                      one_port_operational_cost = '$\mathcal{O}_{n,t}$',
                      co2_cost = '$\mathcal{E}_{n,t}$',
                      nodal_demand_cost = '$\lambda_{n,t} \, d_{n,t}$')


n = ntl.test.get_network_ac_dc()
n.carriers = n.carriers.drop('battery')
# n.generators['p_nom_min'] = 0

tag = ''
relax_co2 = False
if relax_co2:
    n.global_constraints = n.global_constraints.drop('co2_limit')
    tag = '_relaxed_co2'

sns = n.snapshots
n.lopf(pyomo=False, keep_shadowprices=True)
n.carriers.loc['wind', 'color'] = 'steelblue'
n.buses.loc[n.buses.index.str.contains('DC'), 'y'] += 1
n.buses.loc[n.buses.index.str.contains('DC'), 'x'] -= .5

# %% Network plot

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,5))
bus_sizes = n.generators.groupby(['bus', 'carrier']).p_nom_opt.sum()
plot = n.plot(bus_sizes=bus_sizes/1e4, line_widths=n.lines.s_nom_opt/1e2,
       link_widths=n.links.p_nom_opt/1e2, ax=ax)
ntl.plot.annotate_bus_names(n, ax, shift=(0,0.3),
                            transform=ccrs.PlateCarree(),
                            bbox='fancy')
ax.legend(*ntl.plot.handles_labels_for(n.carriers.color), loc='upper left')
fig.canvas.draw(); fig.tight_layout()
fig.savefig(f'figures/network{tag}.png')

# %%
ds = ntl.allocate_flow(n, method='ebe', q=0, aggregated=False)
dc = nodal_demand_cost(n).rename(bus='payer')
ca = ntl.allocate_cost(n, method=ds, q=0)
ca = ca.sum([d for d in ca.dims if d not in ['payer', 'snapshot']])

# %% Allocation plot

color = pd.Series({'one_port_operational_cost': 'darkkhaki',
                   'co2_cost': 'tomato',
                   'one_port_investment_cost': 'palevioletred',
                   'branch_investment_cost': 'mediumaquamarine',
                   'nodal_demand_cost': 'cadetblue'})

fig, axes = plt.subplots(3,3, figsize=(12,7), sharex=True)

for sn, ax in zip(n.snapshots, axes.ravel()):
    ca_df = ca.sel(snapshot=sn, drop=True).to_dataframe()
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
labels = [to_symbol_dict[l] for l in labels]
fig.legend(handles[:-1], labels[:-1], ncol=4, frameon=False, loc='lower center',
           bbox_to_anchor=(0.3, 1), fontsize='large', title='Flow Based (left bars)')
fig.legend(handles[-1:], labels[-1:], ncol=1, frameon=False, loc='lower center',
           bbox_to_anchor=(0.7, 1), fontsize='large', title='LMP Based (right bars)')
fig.tight_layout()
fig.savefig(f'figures/compare_allocation{tag}.png', bbox_inches='tight')

# %% Nodal payments

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(5,4))
bus_sizes = ca.sum('snapshot').to_dataframe().stack()
n.plot(ax=ax, bus_sizes=bus_sizes/bus_sizes.sum()*3, bus_colors=color)
handles, labels = ntl.plot.handles_labels_for(color[:-1])
labels = ['$\sum_t' + to_symbol_dict[l][1:] for l in labels]
ax.legend(handles, labels, loc='upper left')
fig.canvas.draw(); fig.tight_layout()
fig.savefig(f'figures/nodal_payments{tag}.png')


#%% opex flow
from netallocation.cost import (allocate_one_port_operational_cost,
                                allocate_one_port_investment_cost)
from netallocation.plot import pos_neg_buscolors as buscolors
norm = lambda ds: ds/ ds.abs().sum()
t = n.snapshots[0]
F = ds.virtual_flow_pattern
P = ds.virtual_injection_pattern

opex_flow = allocate_one_port_operational_cost(F, n, dim='bus')\
            .sel(snapshot=t).sum(['source_carrier', 'bus']).to_series()

opex_injection = allocate_one_port_operational_cost(P, n, dim='injection_pattern')\
                 .sel(snapshot=t).sum(['source_carrier', 'injection_pattern'])\
                 .to_series()

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.EqualEarth()}, figsize=(8,6))
n.plot(flow=norm(opex_flow) * 100, bus_sizes=norm(opex_injection.abs()) * 0.5,
       bus_colors=buscolors(opex_injection), title=f'OPEX flow Hour {t}',
       ax=ax, line_colors='teal', link_colors='teal')
fig.canvas.draw(); fig.tight_layout()
fig.savefig(f'figures/opex_flow{tag}.png')


#%% capex flow
t = n.snapshots[5]

capex_flow = allocate_one_port_investment_cost(F, n, dim='bus')\
            .sel(snapshot=t).sum(['source_carrier', 'bus']).to_series()

capex_injection = allocate_one_port_investment_cost(P, n, dim='injection_pattern')\
                 .sel(snapshot=t).sum(['source_carrier', 'injection_pattern'])\
                 .to_series()

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.EqualEarth()}, figsize=(8,6))
n.plot(flow=norm(capex_flow) * 100, bus_sizes=norm(capex_injection.abs()) * 0.5,
        bus_colors=buscolors(capex_injection), title=f'CAPEX flow Hour {t}',
        ax=ax, line_colors='teal', link_colors='teal')
fig.canvas.draw(); fig.tight_layout()
fig.savefig(f'figures/capex_flow{tag}.png')

