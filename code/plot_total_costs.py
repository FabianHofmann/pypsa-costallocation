#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 21:56:19 2020

@author: fabian
"""
import pypsa
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from config import color, to_total_symbol
import numpy as np
from pypsa.descriptors import nominal_attrs

plt.rc('font', family='serif')
plt.rc("text", usetex=False)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_total_costs', nname='acdc')

n = pypsa.Network(snakemake.input.network)
ca = xr.open_dataset(snakemake.input.costs)



stacked_cost = ca.sum().to_array().to_series()


if 'lv_limit' in n.global_constraints.index:
    for c in n.branch_components:
        if n.df(c).empty: continue
        mu_upper = (n.global_constraints.at['lv_limit', 'mu'] * n.df(c).length)
        nom = nominal_attrs[c]
        n.df(c)['mu_upper_'+nom] += mu_upper

stacked_cost['subsidy_cost'] = 0
stacked_cost['scarcity_cost'] = 0
for c in ['Generator', 'StorageUnit', 'Line', 'Link']:
    nom = nominal_attrs[c]
    stacked_cost['subsidy_cost'] += n.df(c)[nom + '_min'] @ n.df(c)['mu_lower_' + nom]
    stacked_cost['scarcity_cost'] += (n.df(c)[nom + '_max'].replace(np.inf, 0) @
                                      n.df(c)['mu_upper_' + nom])


assert ((n.objective_constant + n.objective)/stacked_cost.sum()).round(2) == 1

# %%
fig, ax = plt.subplots(figsize=[3,5])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

color = color[stacked_cost.index]
label = pd.Series(to_total_symbol)[stacked_cost.index]

stacked_cost.to_frame('').T.div(1e9).plot.bar(color=color,  stacked=True, ax=ax, zorder=4)
ax.legend(labels=label, loc='lower center', bbox_to_anchor=(0.5, 1), ncol=3, frameon=False)
ax.grid(axis='y', linestyle='dashed', color='gray', zorder=1, alpha=.5)
ax.set_ylabel('Total Cost [Billion €]')
fig.tight_layout()

fig.savefig(snakemake.output[0], bbox_inches='tight')