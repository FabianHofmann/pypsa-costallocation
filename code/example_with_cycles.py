#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 18:44:55 2020

@author: fabian
"""

# Add lp-costallocation  to path,
import netallocation as ntl
import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('text.latex', preamble=r'\usepackage{accents}')


if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('plot_example_with_cycles')


# %% create network

n = pypsa.Network()
n.madd("Bus", range(1,4), x=[0,0,0.5], y=[0,0.7, 0.35])
n.madd('Generator', [1, 3], bus=[1, 3], p_nom=[50, 50],
       marginal_cost=[6, 4], carier='gen')
n.madd("Load", [1,2], bus=[1,2], p_set=[30,50], carrier='load')
n.madd('Line', ['1-2', '3-1', '3-2',],
       bus0=[1, 3, 3,], bus1=[2, 1, 2,],
       s_nom=[60, 30, 30], x=[0.01, 0.01, 0.01], )

n.lopf(pyomo=False, solver_name='gurobi', keep_shadowprices=True,
        keep_references=True)


# %% get quantities

production = n.generators_t.p.loc['now']
demand = n.loads_t.p.loc['now']
lmp = n.buses_t.marginal_price.loc['now']
l_diff = - ntl.cost.locational_market_price_diff(n).squeeze().to_series()['Line']
l_cycle = ntl.cost.cycle_constraint_price(n).squeeze().to_series()['Line']


# %% plot figure
fig, ax = plt.subplots(figsize=(5,6))
n.plot(ax=ax, geomap=False, flow='now', line_widths=0.04, bus_sizes=5e-3)

# bus boxes
bbox = dict(facecolor='lightgrey', alpha=.15, boxstyle='round', pad=0.5)
textkwargs = dict(size=13, va="top", bbox=bbox, zorder=8)
boxes = {bus: r'\hspace{-14pt}' for bus in n.buses.index}

for bus in n.buses.index:

    load_i = n.loads.query('bus == @bus').index
    if not load_i.empty:
        d = n.loads_t.p[load_i.item()]['now']
        boxes[bus] += r'd = %i MW\\'%d

    gen_i = n.generators.query('bus==@bus').index
    if not gen_i.empty:
        g = n.generators_t.p[gen_i.item()]['now']
        boxes[bus] += r'g = %i MW\\'%round(g,1)

    l = n.buses_t.marginal_price[bus].now
    l = '%i'%l if l == round(l,0) else '%.2f'%l
    boxes[bus] += r'$\lambda$ = %s €/MWh \\'%l

    boxes[bus] = boxes[bus][:-2]

ax.text(n.buses.x[0]-.12, n.buses.y[0]+.04, boxes['1'], ha='right', **textkwargs)
ax.text(n.buses.x[1]-.12, n.buses.y[1]+.04, boxes['2'], ha='right', **textkwargs)
ax.text(n.buses.x[2]-.1, n.buses.y[2]-.11, boxes['3'], ha='left', **textkwargs)


# line boxes
bbox = dict(facecolor='indianred', alpha=.15, boxstyle='round', pad=0.5)
textkwargs = dict(size=13, bbox=bbox, zorder=8)
boxes = {bus: r'\hspace{-14pt}' for bus in n.lines.index}

for line in n.lines.index:

    f = n.lines_t.p0[line]['now']
    boxes[line] += r'f = %i MW\\'%f

    l = l_diff[line]
    if l != 0:
        l = '%i'%l if l == round(l,0) else '%.2f'%l
        boxes[line] += r'$\lambda^\mathrm{diff}$ = %s €/MW\\'%l if l != 0 else ''

    l = l_cycle[line]
    if l != 0:
        l = '%i'%l if l == round(l,0) else '%.2f'%l
        boxes[line] += r'$\lambda^\mathrm{KVL}$ = %s €/MW\\'%l if l != 0 else ''

    boxes[line] = boxes[line][:-2]

n.lines['xpos'] = n.lines.bus0.map(n.buses.x)/2 + n.lines.bus1.map(n.buses.x)/2
n.lines['ypos'] = n.lines.bus0.map(n.buses.y)/2 + n.lines.bus1.map(n.buses.y)/2

ax.text(n.lines.xpos[0]-.1, n.lines.ypos[0], boxes['1-2'], ha='right', va='center', **textkwargs)
ax.text(n.lines.xpos[1]+.1, n.lines.ypos[1]+-.1, boxes['3-1'], ha='center', va='top', **textkwargs)
ax.text(n.lines.xpos[2]+.2, n.lines.ypos[2]+.1, boxes['3-2'], ha='center', va='bottom', **textkwargs)

ntl.plot.annotate_bus_names(n, ax, shift=0, color='w', size=18)
fig.tight_layout()
fig.savefig(snakemake.output.network, bbox_inches='tight')

# %% revenue allocation

C = ntl.allocate_revenue(n, with_cycle_prices=True)

to_index = pd.Series(dict(payer='n', receiver_node='s', receiver_branch='$\ell$'))
payoff = C.sum(['snapshot', 'receiver_carrier'])\
           .assign_coords(receiver_branch=['1-2', '1-3', '2-3'])\
           .rename(to_index).transpose(..., 'n') \
           .drop_sel(n='3')
vmax = payoff.to_array().max()

fig, axes = plt.subplots(1, 2, figsize=(4,3), sharey=True,
                         gridspec_kw={'width_ratios':[0.5,.5,]})
for i, v in enumerate(payoff):
    ax = axes[i]
    df = payoff[v].to_pandas().T
    annot = df.round(0).astype(int).astype(str) + ' €'
    annot = annot.where(annot != '0 €', "-")
    sns.heatmap(df, cmap='PRGn', ax=ax, annot=annot, fmt='s',
                cbar=False, vmin=-vmax, vmax=vmax, linewidths=1, square=True,
                annot_kws={'horizontalalignment': 'left'})
    if i:
        ax.set_ylabel('')
    # ax.set_title(v.replace("_", " ").title())
    ax.tick_params(left=False, bottom=False)
fig.tight_layout(w_pad=.4)
fig.savefig(snakemake.output.payoff_kvl, bbox_inches='tight')

# %% revenue allocation wo KVL prices

C = ntl.allocate_revenue(n, with_cycle_prices=False)

to_index = pd.Series(dict(payer='n', receiver_node='s', receiver_branch='$\ell$'))
payoff = C.sum(['snapshot', 'receiver_carrier'])\
           .assign_coords(receiver_branch=['1-2', '1-3', '2-3'])\
           .rename(to_index).transpose(..., 'n') \
           .drop_sel(n='3')
vmax = payoff.to_array().max()

fig, axes = plt.subplots(1, 2, figsize=(4,3), sharey=True,
                         gridspec_kw={'width_ratios':[0.5,.5,]})
for i, v in enumerate(payoff):
    ax = axes[i]
    df = payoff[v].to_pandas().T
    annot = df.round(0).astype(int).astype(str) + ' €'
    annot = annot.where(annot != '0 €', "-")
    sns.heatmap(df, cmap='PRGn', ax=ax, annot=annot, fmt='s',
                cbar=False, vmin=-vmax, vmax=vmax, linewidths=1, square=True,
                annot_kws={'horizontalalignment': 'left'})
    if i:
        ax.set_ylabel('')
    # ax.set_title(v.replace("_", " ").title())
    ax.tick_params(left=False, bottom=False)
fig.tight_layout(w_pad=.4)
fig.savefig(snakemake.output.payoff, bbox_inches='tight')

