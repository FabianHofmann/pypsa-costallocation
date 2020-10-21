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
        snakemake = mock_snakemake('plot_total_cost_comparison', nname_wo_field='de10')

ngf = pypsa.Network(snakemake.input.network_gf)
nbf = pypsa.Network(snakemake.input.network_bf)
cagf = xr.open_dataset(snakemake.input.costs_gf)
cabf = xr.open_dataset(snakemake.input.costs_bf)

# %%

TC_gf = ngf.objective_constant + ngf.objective
TC_bf = nbf.objective_constant + nbf.objective

stacked_cost = pd.concat([cagf.sum().to_array().to_series(),
                          cabf.sum().to_array().to_series()],
                         keys=['Greenfield', 'Brownfield'], axis=1).T

stacked_cost['remaining_cost'] =  [TC_gf, TC_bf] - stacked_cost.sum(1)

# %%
fig, ax = plt.subplots(figsize=[5,5])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

color = color[stacked_cost.columns]
label = pd.Series(to_total_symbol)[stacked_cost.columns]

stacked_cost.div(1e9).plot.bar(color=color,  stacked=True, ax=ax, zorder=4)
ax.legend(labels=label, loc='lower center', bbox_to_anchor=(0.5, 1), ncol=5, frameon=False)
ax.grid(axis='y', linestyle='dashed', color='gray', zorder=1, alpha=.5)
ax.set_ylabel('Total Revenue [Billion €]')

fig.savefig(snakemake.output[0])