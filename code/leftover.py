#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 19:41:23 2020

@author: fabian
"""


from netallocation.cost import locational_market_price
from netallocation.grid import (PTDF, Incidence, power_demand, network_injection,
                                power_production)
from netallocation.linalg import dedup_axis, dot
from numpy import eye
import xarray as xr

norm = lambda ds: ds/ds.sum()

m = 'DE0 37'
t = n.snapshots[1]
sel = dict(bus=m, snapshot=t)

k = norm(power_production(n, t))
# k = xr.zeros_like(k) + norm(np.random.rand(*k.shape))
H = PTDF(n) - PTDF(n) @ k
K = Incidence(n)
y = locational_market_price(n, t)
d = ntl.power_demand(n, t).sel(bus=m)

line_mu_up = n.lines_t.mu_upper.loc[t]
line_mu_lo = n.lines_t.mu_lower.loc[t]
lhs = y.sel(bus=m)
rhs = -(line_mu_up - line_mu_lo) @ H.sel(bus=m) + y @ k
print(lhs.item(), rhs.item())

# test equations
o = n.generators.marginal_cost[s] + ntl.cost.nodal_co2_price(n, t).loc[m, carrier]
mu_up = n.generators_t.mu_upper.loc[t, s]
mu_lo = n.generators_t.mu_lower.loc[t, s]
line_mu_up = n.lines_t.mu_upper.loc[t]
line_mu_lo = n.lines_t.mu_lower.loc[t]
rhs  = ((line_mu_up - line_mu_lo) @ H.sel(bus=m))
print(lhs.item(), rhs.item())

lhs = - y.sel(bus=m) * K.sel(bus=m) @ H.sel(bus=m)
wo_m = n.buses.drop(m).index
rhs = ((line_mu_up - line_mu_lo) * H.sel(bus=m)).sum() +\
       y.sel(bus=wo_m) @ K.sel(bus=wo_m) @ H.sel(bus=m)
print(lhs.item(), rhs.item())


# %% test standard equation

from netallocation.cost import reindex_by_bus_carrier, nodal_co2_price
c = 'Generator'
(
 - reindex_by_bus_carrier(n.df(c).marginal_cost, c, n) \
 - nodal_co2_price(n, t) \
 - reindex_by_bus_carrier(n.pnl(c).mu_upper.loc[t], c, n) \
 + reindex_by_bus_carrier(n.pnl(c).mu_lower.loc[t], c, n) \
 + locational_market_price(n, t)
 )

#%% comparison plot bars

fig, axes = plt.subplots(3, 1, figsize=(12, 7), sharex=True)

for sn, ax in zip(n.snapshots, axes.ravel()):
    ca_df = payments.sel(snapshot=sn, drop=True).to_dataframe()
    dc_df = dc.sel(snapshot=sn, drop=True).to_dataframe()
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
fig.savefig(f"../figures/compare_allocation{tag}.png", bbox_inches="tight")


