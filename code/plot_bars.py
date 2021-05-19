#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 10:12:28 2020

@author: fabian
"""

import xarray as xr
import matplotlib.pyplot as plt
import pypsa
from config import color, to_symbol, sink_dims, source_carrier
import netallocation as ntl


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_bars', nname='test-de10')


n = pypsa.Network(snakemake.input.network)
payments = xr.open_dataset(snakemake.input.payments).sum(source_carrier)

if snakemake.config['alloc_to_load_only']:
    demand = ntl.power_demand(n, per_carrier=True)\
                .sel(carrier='Load', drop=True).rename(bus='sink')
else:
    demand = ntl.power_demand(n).rename(bus='sink')

dc = demand * ntl.cost.locational_market_price(n).rename(bus='sink')

fig, axes = plt.subplots(3, 1, figsize=(12, 7), sharex=True)

for sn, ax in zip(n.snapshots, axes.ravel()):
    for spine in ['right', 'top', 'bottom']:
        ax.spines[spine].set_visible(False)
    ca_df = payments.sel(snapshot=sn, drop=True).to_dataframe()
    dc_df = dc.sel(snapshot=sn, drop=True).to_dataframe(name='nodal_demand_cost')
    d = dict(stacked=True, zorder=3, width=0.3, legend=False, ax=ax)

    ca_df.plot.bar(position=1.1, color=color[ca_df.columns], **d)
    dc_df.plot.bar(position=0, color=color[dc_df.columns], **d)

    ax.set_xlim(left=-0.5)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(5, 2))
    ax.set_xlabel("")
    ax.set_title("$t = %i$" % (sn.hour + 1))
    ax.grid(axis="y", zorder=1, linestyle="dashed")
    ax.hlines([0], *ax.get_xlim(),  ls='dashed', alpha=0.7, color='grey', lw=1)

handles, labels = ax.get_legend_handles_labels()
labels = [to_symbol[l] for l in labels]
fig.legend(
    handles[:-1],
    labels[:-1],
    ncol=5,
    frameon=False,
    loc="lower center",
    bbox_to_anchor=(0.3, 1),
    fontsize="large",
    title="Cost Allocation (left bars)",
)
fig.legend(
    handles[-1:],
    labels[-1:],
    ncol=1,
    frameon=False,
    loc="lower center",
    bbox_to_anchor=(0.7, 1),
    fontsize="large",
    title="Nodal Expenditures (right bars)",
)
fig.tight_layout()
fig.savefig(snakemake.output[0], bbox_inches="tight")


