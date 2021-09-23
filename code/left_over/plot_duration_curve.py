#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 10:34:01 2020

@author: fabian
"""
import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import config


if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('plot_capex_duration_curve', nname='de50')


n = pypsa.Network(snakemake.input.network)


c = 'Generator'
gen = (n.pnl(c).p * n.pnl(c).mu_upper).groupby(n.df(c).carrier, axis=1).sum()

c = 'StorageUnit'
mu = (n.pnl(c).mu_state_of_charge / n.df(c).efficiency_dispatch
          + n.pnl(c).mu_upper_p_dispatch + n.pnl(c).mu_lower_p_dispatch)
sto = (n.pnl(c).p_dispatch * mu).groupby(n.df(c).carrier, axis=1).sum().clip(lower=0)

capex = pd.concat([gen, sto], axis=1)
# %%
fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (9, 3), gridspec_kw={'width_ratios':[1, 5]},
                               sharey=True)

duration = (capex.reindex(index=capex.sum(1).sort_values(ascending=False).index,
                          columns=capex.sum().sort_values(ascending=False).index)
            .reset_index(drop=True))

# mimic a bar plot with area being filled close to the next number
duration = pd.concat([duration, duration.rename(lambda x: x+.99)]).sort_index()

# split the x axis here
xsplit = 10


duration[:xsplit*2].plot.area(stacked=True, ax=ax1, legend=False,
                                color=n.carriers.color[duration.columns])

duration[xsplit*2:].plot.area(stacked=True, ax=ax2, legend='reverse',
                                color=n.carriers.color[duration.columns])

handles = ax2.get_legend_handles_labels()[0]
ax2.legend(handles[::-1], n.carriers.nice_name[duration.columns][::-1],
           ncol=3, frameon=False)


ax1.set_ylabel('Allocated CAPEX [€]')
ax2.set_xlabel('time [h]')
ax1.set_yscale('log')
ax1.set_ylim(bottom=10)
ax1.set_xlim(right=xsplit)
ax2.set_xlim(left=xsplit)

locs, labels = plt.xticks()

ticks = ax2.get_xticks()
ticks[0] = xsplit
ax2.set_xticks(ticks)


ax1.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

ax2.get_yaxis().set_visible(False)

ax1.spines['right'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)


fig.tight_layout()
fig.savefig(snakemake.output[0], bbox_inches='tight')
