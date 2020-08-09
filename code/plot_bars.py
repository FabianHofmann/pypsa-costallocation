#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 10:12:28 2020

@author: fabian
"""

import xarray as xr
import matplotlib.pyplot as plt
import pypsa
from config import color, to_symbol, sink_dims
import netallocation as ntl


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_bars', nname='acdc',
                                   method='ptpf', power='net')

load_only = True

n = pypsa.Network(snakemake.input.network)
payments = xr.open_dataset(snakemake.input.payments)
demand = ntl.power_demand(n, per_carrier=True).rename(sink_dims)
dc = demand * ntl.cost.locational_market_price(n).rename(bus='sink')

if load_only:
    dc = dc.sel(sink_carrier='Load', drop=True)
    payments = payments.sel(sink_carrier='Load', drop=True)
else:
    dc = dc.sum('sink_carrier')
    payments = payments.sum('sink_carrier')

fig, axes = plt.subplots(3, 1, figsize=(12, 7), sharex=True)

for sn, ax in zip(n.snapshots, axes.ravel()):
    ca_df = payments.sel(snapshot=sn, drop=True).to_dataframe()
    dc_df = dc.sel(snapshot=sn, drop=True).to_dataframe(name='nodal_demand_cost')
    d = dict(stacked=True, zorder=3, width=0.3, legend=False, ax=ax)

    ca_df.plot.bar(position=1.1, color=color[ca_df.columns], **d)
    dc_df.plot.bar(position=0, color=color[dc_df.columns], **d)

    ax.set_xlim(left=-0.5)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(5, 2))
    ax.set_xlabel("")
    ax.set_title(f"$t = %i$" % (sn.hour + 1))
#     ax.grid(axis="y", zorder=1, linestyle="dashed")

handles, labels = ax.get_legend_handles_labels()
labels = [to_symbol[l] for l in labels]
fig.legend(
    handles[:-1],
    labels[:-1],
    ncol=4,
    frameon=False,
    loc="lower center",
    bbox_to_anchor=(0.3, 1),
    fontsize="large",
    title="Flow Based (left bars)",
)
fig.legend(
    handles[-1:],
    labels[-1:],
    ncol=1,
    frameon=False,
    loc="lower center",
    bbox_to_anchor=(0.7, 1),
    fontsize="large",
    title="LMP Based (right bars)",
)
fig.tight_layout()
fig.savefig(snakemake.output[0], bbox_inches="tight")


