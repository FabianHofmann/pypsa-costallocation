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
        snakemake = mock_snakemake('plot_total_costs', nname_wo_field='de10')

ngf = pypsa.Network(snakemake.input.network_gf)
nbf = pypsa.Network(snakemake.input.network_bf)
cagf = xr.open_dataset(snakemake.input.costs_gf)
cabf = xr.open_dataset(snakemake.input.costs_bf)

# %%
fig, ax = plt.subplots()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

TC_gf = ngf.objective_constant + ngf.objective
TC_bf = nbf.objective_constant + nbf.objective
TC = [TC_gf, TC_bf]

stacked_gf = (cagf.sum().to_array().to_series().reindex_like(color).dropna()
              .rename(to_total_symbol).rename('').to_frame().T)
stacked_bf = (cabf.sum().to_array().to_series().reindex_like(color).dropna()
              .rename(to_total_symbol).rename('').to_frame().T)
stacked_cost = pd.concat([stacked_gf, stacked_bf], keys=['Greenfield', 'Brownfield'])

R = TC - stacked_cost.sum(1)
stacked_cost[r'$\mathcal{R}$'] = R

stacked_cost.plot.bar(color=color, stacked=True, ax=ax, zorder=4)
ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1), ncol=5, frameon=False)
ax.grid(axis='y', linestyle='dashed', color='gray', zorder=1, alpha=.5)

fig.savefig(snakemake.output[0])