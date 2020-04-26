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


to_symbol_dict = dict(one_port_investment_cost = '$\mathcal{C}^{G}_{m,t}$',
                      branch_investment_cost = '$\mathcal{C}^{F}_{m,t}$',
                      one_port_operational_cost = '$\mathcal{O}_{m,t}$',
                      co2_cost = '$\mathcal{E}_{m,t}$',
                      nodal_demand_cost = '$\sum_a \, \lambda_{m,t} \, d_{m,a,t}$')


n = pypsa.Network('elec_s_50_ec_lvopt_Co2L-3H.nc')
n.calculate_dependent_values()
n.snapshot_weightings[:] = 1
n.generators['p_nom_extendable'] = True
n.lines['s_nom_min'] = 0
n.generators['p_nom_max'] = np.inf
n.mremove('StorageUnit', n.storage_units.index)
n.mremove('Link', n.links.index)
n.mremove('Line', n.lines.query('x_pu == inf').index)
n.remove('Bus', 'DE0 0')

T = 10
n.set_snapshots(n.snapshots[:T])
n.global_constraints.loc['CO2Limit', 'constant'] *= T/8700
n.lopf(pyomo=False, keep_shadowprices=True, solver_name='gurobi')
tag = '_large'

# %% reproduce analytic

from netallocation.cost import locational_market_price
from netallocation.grid import PTDF, Incidence, power_demand, network_injection
from netallocation.linalg import dedup_axis, dot
from numpy import eye

norm = lambda ds: ds/ds.sum()

m = 'DE0 9'
carrier = 'solar'
s = m + ' ' + carrier
t = n.snapshots[2]

slack = norm(network_injection(n, t).clip(max=0))
k = PTDF(n) @ slack
H = PTDF(n) - k
K = Incidence(n)
y = locational_market_price(n, t)
lhs = - (y @ K @ H.sel(bus=m))

o = n.generators.marginal_cost[s] + ntl.cost.nodal_co2_price(n, t).loc[m, carrier]
mu_up = n.generators_t.mu_upper.loc[t, s]
mu_lo = n.generators_t.mu_lower.loc[t, s]
line_mu_up = n.lines_t.mu_upper.loc[t]
line_mu_lo = n.lines_t.mu_lower.loc[t]
rhs  = ((line_mu_up - line_mu_lo) @ H.sel(bus=m))

print(lhs.item(), rhs.item())

lhs = - y.sel(bus=m) * K.sel(bus=m) @ H.sel(bus=m)
wo_m = n.buses.drop(m).index
rhs = ((line_mu_up - line_mu_lo) * H.sel(bus=m)).sum() +\
       y.sel(bus=wo_m) @ K.sel(bus=wo_m) @ H.sel(bus=m)
print(lhs.item(), rhs.item())


N = len(n.buses)
D = dedup_axis(dot(K, H), ('bus_0', 'bus_1'))
print(D.to_pandas().round(3))
# %%
p = D.sel(bus_1=m).to_series()
f = H.sel(bus=m).to_series()
c = pd.Series('red', p.index).where(p<0, 'blue')
n.plot(flow=f*50, bus_sizes=p.abs(), bus_colors=c, bus_alpha=0.5); ntl.plot.annotate_bus_names(n, shift=0)

# %% test standard equation

from netallocation.cost import reindex_by_bus_carrier, nodal_co2_price
c = 'Generator'
(
 - reindex_by_bus_carrier(n.df(c).marginal_cost, c, n) \
 - nodal_co2_price(n, t) \
 - reindex_by_bus_carrier(n.pnl(c).mu_upper.loc[t], c, n) \
 + reindex_by_bus_carrier(n.pnl(c).mu_lower.loc[t], c, n) \
 + locational_market_price(n, t)
 )


# %% Network plot

# fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(6,5))
# bus_sizes = n.generators.groupby(['bus', 'carrier']).p_nom_opt.sum()
# plot = n.plot(bus_sizes=bus_sizes/1e5, line_widths=n.lines.s_nom_opt/1e3,
#        link_widths=n.links.p_nom_opt/1e3, ax=ax)
# ax.legend(*ntl.plot.handles_labels_for(n.carriers.color), loc='center left',
#           bbox_to_anchor=(1, 0.5), ncol=1)
# fig.canvas.draw(); fig.tight_layout()
# fig.savefig(f'figures/network{tag}.png')

# # %%

# ds = ntl.allocate_flow(n, method='ebe', q=0, direct=False)
# pr = nodal_production_revenue(n).rename(bus='payer')
# dc = nodal_demand_cost(n).rename(bus='payer')
# ca = ntl.allocate_cost(n, method=ds, q=0)
# ca = ca.sum([d for d in ca.dims if d not in ['payer', 'snapshot']])
# dc = nodal_demand_cost(n).rename(bus='payer')

# # %% Allocation plot

# color = pd.Series({'one_port_operational_cost': 'darkkhaki',
#                    'co2_cost': 'tomato',
#                    'one_port_investment_cost': 'palevioletred',
#                    'branch_investment_cost': 'mediumaquamarine',
#                    'nodal_demand_cost': 'cadetblue'})

# fig, axes = plt.subplots(3,1, figsize=(12,7), sharex=True)

# for sn, ax in zip(n.snapshots, axes.ravel()):
#     ca_df = ca.sel(snapshot=sn, drop=True).to_dataframe()
#     dc_df = dc.sel(snapshot=sn, drop=True).to_dataframe()
#     d = dict(stacked=True, zorder=3, width=.3, legend=False, ax=ax)

#     ca_df.plot.bar(position=1.1, color=color[ca_df.columns], **d)
#     dc_df.plot.bar(position=0, color=color[dc_df.columns], **d)

#     ax.set_xlim(left=-.5)
#     ax.ticklabel_format(axis="y", style="sci", scilimits=(5,2))
#     ax.set_xlabel('')
#     ax.set_title(f'$t = %i$'%(sn.hour+1))
#     ax.grid(axis='y', zorder=1, linestyle='dashed')

# handles, labels = ax.get_legend_handles_labels()
# labels = [to_symbol_dict[l] for l in labels]
# fig.legend(handles[:-1], labels[:-1], ncol=4, frameon=False, loc='lower center',
#            bbox_to_anchor=(0.3, 1), fontsize='large', title='Flow Based (left bars)')
# fig.legend(handles[-1:], labels[-1:], ncol=1, frameon=False, loc='lower center',
#            bbox_to_anchor=(0.7, 1), fontsize='large', title='LMP Based (right bars)')
# fig.tight_layout()
# fig.savefig(f'figures/compare_allocation{tag}.png', bbox_inches='tight')

# #%% Nodal payments
# fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(5,4))
# bus_sizes = ca.sum('snapshot').to_dataframe().stack()
# n.plot(ax=ax, bus_sizes=bus_sizes/bus_sizes.sum()*3, bus_colors=color)
# handles, labels = ntl.plot.handles_labels_for(color[:-1])
# labels = ['$\sum_t' + to_symbol_dict[l][1:] for l in labels]
# ax.legend(handles, labels, loc='upper left')
# fig.canvas.draw(); fig.tight_layout()
# fig.savefig(f'figures/nodal_payments{tag}.png')


# #%% opex flow
# norm = lambda ds: ds/ ds.abs().sum()
# t = 0
# opex = ntl.cost.allocate_one_port_operational_cost(ds, n).sel(snapshot=sns[t])
# opex = opex.sum(['source_carrier'])
# opex_flow = opex.peer_on_branch_to_peer.sum(['source', 'sink']).to_series()
# opex_injection = opex.peer_to_peer.sum('sink').rename(source='bus') -\
#                  opex.peer_to_peer.sum('source').rename(sink='bus')
# opex_injection = opex_injection.to_series()
# buscolors = pd.Series('indianred', n.buses.index).where(opex_injection<0, 'steelblue')

# fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8,6))
# n.plot(flow=norm(opex_flow) * 100, bus_sizes=norm(opex_injection.abs()), bus_colors=buscolors,
#        title=f'OPEX flow Hour {t}', ax=ax, line_colors='teal', link_colors='teal')
# fig.canvas.draw(); fig.tight_layout()
# fig.savefig(f'figures/opex_flow{tag}.png')


# #%% capex flow
# t = 5
# capex = ntl.cost.allocate_one_port_investment_cost(ds, n).sel(snapshot=sns[t])
# capex = capex.sum(['source_carrier'])
# capex_flow = capex.peer_on_branch_to_peer.sum(['source', 'sink']).to_series()
# capex_injection = capex.peer_to_peer.sum('sink').rename(source='bus') -\
#                   capex.peer_to_peer.sum('source').rename(sink='bus')
# capex_injection = capex_injection.to_series()
# buscolors = pd.Series('indianred', n.buses.index).where(capex_injection<0, 'steelblue')
# fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8,6))
# n.plot(flow=norm(capex_flow) * 100, bus_sizes=norm(capex_injection.abs()),
#        bus_colors=buscolors,
#        title=f'CAPEX flow Hour {t}', ax=ax, line_colors='teal', link_colors='teal')
# fig.canvas.draw(); fig.tight_layout()
# fig.savefig(f'figures/capex_flow{tag}.png')

