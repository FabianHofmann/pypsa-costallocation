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

plt.rc('font', family='serif')
plt.rc("text", usetex=False)


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_total_costs', nname='de10bf')

n = pypsa.Network(snakemake.input.network)
ca = xr.open_dataset(snakemake.input.costs)

TC = n.objective_constant + n.objective

stacked_cost = ca.sum().to_array().to_series()
stacked_cost['remaining_cost'] =  TC - stacked_cost.sum()

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