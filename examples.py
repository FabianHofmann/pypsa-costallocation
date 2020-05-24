#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:04:53 2020

@author: fabian
"""

import netallocation as ntl
from netallocation.plot import handles_labels_for
import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from config import to_symbol_dict, color
plt.rc('text', usetex=False)

tag = ''
add_gen_args = {}

# tag = '_constraint_capacity'
# add_gen_args = {'p_nom_max': 80}


n = pypsa.Network()
n0 = '0'
n1 = '1'
d0 = 50
d1 = 80
o0 = 50
o1 = 200
c0 = 500
c1 = 500
c_l = 100

n.madd('Bus', [n0, n1], x=[0, 1], y=[0, 0.5])
n.madd('Load', [0, 1], bus=[n0, n1], p_set=[d0, d1])
n.madd('Generator', [0, 1], bus=[n0, n1], p_nom_extendable=True,
       marginal_cost=[o0, o1], capital_cost=[c0, c1], **add_gen_args)
n.madd('Line', [0], s_nom_extendable=True, x=0.01, bus0=[n0], bus1=[n1],
       capital_cost=c_l)
n.lopf(pyomo=False, keep_shadowprices=True)

sn = 'now'

G0 = int(n.generators.p_nom_opt[0])
G1 = int(n.generators.p_nom_opt[1])
F = int(n.lines.s_nom_opt[0])

g0 = int(n.generators_t.p['0']['now'])
g1 = int(n.generators_t.p['1']['now'])
f = int(n.lines_t.p0['0']['now'])

ds = ntl.allocate_flow(n, method='ebe', q=0, aggregated=False)
dc = ntl.cost.nodal_demand_cost(n).rename(bus='payer')
ca = ntl.allocate_cost(n, method=ds, q=0)
ca_sum = ca.sum([d for d in ca.dims if d not in ['payer', 'snapshot']])
ca_df = ca_sum.sel(snapshot=sn, drop=True).to_dataframe()
dc_df = dc.sel(snapshot=sn, drop=True).to_dataframe()

# %%
fig, ax = plt.subplots()

n.plot(flow=ntl.network_flow(n, sn).to_series()/50, margin=0.3, ax=ax,
       geomap=False, bus_sizes = ca_df.stack().clip(lower=0)/5e6, bus_colors=color)
bbox = dict(facecolor='w', alpha=0.3, edgecolor='grey', boxstyle='round')
textkwargs = dict(size=12, color='darkslategray', ha="left", va="top",
                  bbox=bbox, zorder=8)
bus0_fix = (f'Fixed:  \n\no = {o0} €/MW\nc = {c0} €/MW\nd = {d0} MW')
bus0_opt = (f'Optimized:  \n\ng = {g0} MW\nG = {G0} MW')
bus0_pay = (f'Fix quantities:  \n\no = {o0} €/MW\nc = {c0} €/MW\nd = {d0} MW')
ax.text(n.buses.x[0]-.22, n.buses.y[0]+.26, bus0_fix, **textkwargs)
ax.text(n.buses.x[0], n.buses.y[0]+.26, bus0_opt, **textkwargs)

bus1_fix = (f'Fixed:  \n\no = {o1} €/MW\nc = {c1} €/MW\nd = {d1} MW')
bus1_opt = (f'Optimized:  \n\ng = {g1} MW\nG = {G1} MW')
bus0_pay = (f'Fix quantities:  \n\no = {o0} €/MW\nc = {c0} €/MW\nd = {d0} MW')
ax.text(n.buses.x[1]-.18, n.buses.y[1]-.12, bus1_fix, **textkwargs)
ax.text(n.buses.x[1]+.04, n.buses.y[1]-.12, bus1_opt, **textkwargs)

line_fix = (f'Fixed:  \n\nc = {c_l} €/MW')
line_opt = (f'Optimized:  \n\nf = {F} MW\nF = {F} MW')
bus0_pay = (f'Fix quantities:  \n\no = {o0} €/MW\nc = {c0} €/MW\nd = {d0} MW')
ax.text(n.buses.x.mean()-.2, n.buses.y.mean()+.25, line_fix, **textkwargs)
ax.text(n.buses.x.mean()+.02, n.buses.y.mean()+.25, line_opt, **textkwargs)

ax.legend(*handles_labels_for(color[ca_df.columns].rename(to_symbol_dict)),
          ncol=1, loc='upper left')
ntl.plot.annotate_bus_names(n, ax, shift=0, bbox='fancy')
# fig.tight_layout()
fig.savefig(f'figures/example_network{tag}.png', bbox_inches='tight')

# %%
fig, ax = plt.subplots(figsize=(6,6))
d = dict(stacked=True, zorder=3, width=.3, legend=False, ax=ax)

ca_df.plot.bar(position=1.1, color=color[ca_df.columns], **d)
dc_df.plot.bar(position=0, color=color[dc_df.columns], **d)

ax.set_xlim(left=-.5)
ax.ticklabel_format(axis="y", style="sci", scilimits=(5,2))
ax.set_xlabel('')
ax.grid(axis='y', zorder=1, linestyle='dashed')

handles, labels = ax.get_legend_handles_labels()
labels = [to_symbol_dict[l] for l in labels]
fig.legend(handles[:-1], labels[:-1], ncol=4, frameon=False, loc='lower center',
           bbox_to_anchor=(0.3, 1), fontsize='large', title='Flow Based (left bars)')
fig.legend(handles[-1:], labels[-1:], ncol=1, frameon=False, loc='lower center',
           bbox_to_anchor=(0.7, 1), fontsize='large', title='LMP Based (right bars)')
fig.tight_layout()
fig.savefig(f'figures/example_cost_shares{tag}.png', bbox_inches='tight')


# %% payoff matrix

to_indices = pd.Series(dict(payer='n', receiver_nodal_cost='m',
                            receiver_transmission_cost='$\ell$'))
fig, axes = plt.subplots(1, 3, figsize=(10,3))
L = len(ca.receiver_transmission_cost)
payoff = ca.sum(['snapshot', 'receiver_carrier'])\
           .assign_coords(receiver_transmission_cost=range(L))
payoff = payoff.rename(to_indices)
for i, v in enumerate(payoff):
    ax = axes[i]
    df = payoff[v].to_pandas().T
    sns.heatmap(df, cmap='Blues', ax=ax, annot=True, fmt='10.0f', cbar=False)
    ax.set_title(to_symbol_dict[v])
fig.tight_layout()
fig.savefig(f'figures/example_payoff{tag}.png', bbox_inches='tight')

