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
from netallocation.plot import injection_plot_kwargs
from helpers import adjust_shadowprice
from config import to_static_symbol, color, to_total_symbol

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

path = '../figures/example-old'
net = True
power = 'net' if net else 'gross'
method = 'ptpf'


n = pypsa.Network()
n0 = '1'
n1 = '2'
d0 = 60
d1 = 90
o0 = 50
o1 = 200
c0 = 500
c1 = 500
c_l = 100

n.madd('Bus', [n0, n1], x=[0, 1], y=[0, 0.4])
n.madd('Load', [0, 1], bus=[n0, n1], p_set=[d0, d1])
n.madd('Generator', [0, 1], bus=[n0, n1], p_nom_extendable=True,
       marginal_cost=[o0, o1], capital_cost=[c0, c1], p_nom_max=100)
n.madd('Line', ['1'], s_nom_extendable=True, x=0.01, bus0=[n0], bus1=[n1],
       capital_cost=c_l)
n.lopf(pyomo=False, keep_shadowprices=True, solver_name='gurobi')

sn = 'now'

G0, G1 = n.generators.p_nom_opt
F = n.lines.s_nom_opt.item()
Gupper = n.generators.p_nom_max.unique().item()

g0, g1 = n.generators_t.p.loc['now']
mug0, mug1 = n.generators_t.mu_upper.loc['now']
f = n.lines_t.p0.loc['now'].item()
mul = n.lines_t.mu_upper.loc['now'].item()
l0, l1 = n.buses_t.marginal_price.loc['now']


ds = ntl.allocate_flow(n, method=method, aggregated=net)
ds = ntl.convert.peer_to_peer(ds, n)
ds = ntl.convert.virtual_patterns(ds, n, q=0)

ca = ntl.allocate_cost(n, method=ds)
ca = ca.rename(one_port_operational_cost='generator_operational_cost',
          one_port_investment_cost = 'generator_investment_cost')
ca_sum = ca.sum([d for d in ca.dims if d not in ['payer', 'snapshot']])
C = ca_sum.sel(snapshot=sn, drop=True).to_dataframe()

# %%
fig, ax = plt.subplots(figsize=(12,4))

n.plot(flow=ntl.network_flow(n, sn).to_series()/50, margin=0.3, ax=ax,
       geomap=False, bus_sizes = C.stack().clip(lower=0)/5e6, bus_colors=color)
bbox = dict(facecolor='w', alpha=.15, boxstyle='round', pad=0.5)
textkwargs = dict(size=13, va="top", bbox=bbox, zorder=8)
wpad = 0.04

# first bus box
bus0_fix = (r'Fixed:\vspace{5pt} \\ o = %i €/MWh\\ c = %i €/MW\\ d = %i MW'
            r'\\ $\bar{\textrm{G}}$ = %i MW')%(o0, c0, d0, Gupper)
bus0_opt = (r'Optimized:\vspace{5pt} \\ g = %i MW\\ G = %i MW\\ $\lambda$ = %i €'
            r'\\ $\bar{\mu}$ = %i €')%(g0, G0, l0, mug0)

dy = .3
ax.text(n.buses.x[0]-wpad, n.buses.y[0]+dy, bus0_fix, ha='right', **textkwargs)
ax.text(n.buses.x[0], n.buses.y[0]+dy, bus0_opt, ha='left', **textkwargs)

# second bus box
bus1_fix = (r'Fixed:\vspace{5pt} \\ o = %i €/MWh\\ c = %i €/MW\\ d = %i MW'
            r'\\ $\bar{\textrm{G}}$ = %i MW')%(o1, c1, d1, Gupper)
bus1_opt = (r'Optimized:\vspace{5pt} \\ g = %i MW\\ G = %i MW\\ $\lambda$ = %i €'
            r'\\ $\bar{\mu}$ = %i €')%(g1, G1, l1, mug1)

dy = - .15
dx = .1
ax.text(n.buses.x[1]+dx, n.buses.y[1]+dy, bus1_fix, ha='right', **textkwargs)
ax.text(n.buses.x[1]+dx+wpad, n.buses.y[1]+dy, bus1_opt, ha='left', **textkwargs)

# line box
line_fix = r'Fixed:\vspace{5pt} \\ c = %i €/MW'%(c_l)
line_opt = (r'Optimized:\vspace{5pt} \\ f = %i MW\\ F = %i MW'
            r'\\ $\bar{\mu}$ = %i €')%(f,F,mul)

dy = 0.25
dx = 0.05
ax.text(n.buses.x.mean()+dx-wpad, n.buses.y.mean()+dy, line_fix, ha='right',
        **textkwargs)
ax.text(n.buses.x.mean()+dx, n.buses.y.mean()+dy, line_opt, ha='left',
        **textkwargs)

ax.legend(*handles_labels_for(color[C.columns].rename(to_static_symbol)),
          ncol=3, loc='lower right', frameon=True, fontsize='large')
ntl.plot.annotate_bus_names(n, ax, shift=0, bbox='fancy')
fig.tight_layout()
fig.savefig(f'{path}/example_network.png', bbox_inches='tight')


# %% payoff matrix

to_indices = pd.Series(dict(payer='n', receiver_nodal_cost='s',
                            receiver_transmission_cost='$\ell$'))
fig, axes = plt.subplots(1, 3, figsize=(4.5,2.5), sharey=True,
                          gridspec_kw={'width_ratios':[0.4,.4,.2]})
payoff = ca.sum(['snapshot', 'receiver_carrier'])\
           .assign_coords(receiver_transmission_cost=['1'])\
           .rename(to_indices).transpose(..., 'n')
vmax = payoff.to_array().max()
for i, v in enumerate(payoff):
    ax = axes[i]
    df = payoff[v].to_pandas().T
    annot = df.round(0).div(1e3).astype(int).astype(str) + 'k €'
    sns.heatmap(df, cmap='PRGn', ax=ax, annot=annot, fmt='s',
                cbar=False, vmin=-vmax, vmax=vmax, linewidths=1, square=True,
                annot_kws={'horizontalalignment': 'left'})
    if i:
        ax.set_ylabel('')
    ax.set_title(to_total_symbol[v])
    ax.tick_params(left=False, bottom=False)
fig.tight_layout(w_pad=.4)
fig.savefig(f'{path}/example_payoff_{power}_{method}.png', bbox_inches='tight')

# %% sub flows

# fig, axes = plt.subplots(1, 1, figsize=(5, 5))
bcolor = 'cadetblue'
fcolor = 'darksalmon'

# for bus, ax in zip(n.buses.index, axes):
for bus in n.buses.index:

    fig, ax = plt.subplots(1, 1, figsize=(5, 2.5))
    subflow = ds.virtual_flow_pattern.sel(snapshot='now', bus=bus).to_series()
    p2p = ds.peer_to_peer.sel(sink=bus, snapshot='now').to_series()
    if subflow.item() == 0:
        n.plot(bus_sizes=p2p*3e-4, ax=ax, geomap=False, line_widths=1,
               bus_colors=bcolor, line_colors=fcolor)
    else:
        n.plot(flow=subflow*5e-2, bus_sizes=p2p*3e-4, ax=ax, geomap=False,
               bus_colors=bcolor, line_colors=fcolor)

    bbox2 = dict(facecolor='lightgrey', alpha=1, edgecolor='None', boxstyle='circle')
    ntl.plot.annotate_bus_names(n, ax, shift=0, bbox=bbox2,)

    for b in n.buses.index:
        bbox.update({'facecolor': bcolor, 'edgecolor': 'None', 'pad':.5, 'alpha':1})
        A = r'$A_{%s \rightarrow %s} = %i$ MW'%(b, bus, p2p[b])
        ax.text(*n.buses.loc[b, ['x', 'y']] + [0, 0.2], A, zorder=8, color='white',
                ha='center', va='bottom', bbox=bbox)

    source = p2p.drop(bus).index.item()

    bbox.update({'facecolor': fcolor})
    A = r'$A_{1,%s} = %+i$ MW'%(bus, round(subflow.item()))
    A = A.replace('+', '')
    ax.text(*n.buses[['x', 'y']].mean() + [0, 0.2], A, zorder=8,
            ha='center', va='center', color='white',
            bbox=bbox)
    ax.autoscale()
    fig.tight_layout(h_pad=0, w_pad=0)
    fig.savefig(f'{path}/example_allocation_bus{bus}_{power}_{method}.png',
                bbox_inches='tight')




