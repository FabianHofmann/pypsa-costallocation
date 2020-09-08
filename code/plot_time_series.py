#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 10:34:01 2020

@author: fabian
"""
import pypsa
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('plot_capex_timeseries', nname='acdc')


n = pypsa.Network(snakemake.input.network)

months = mdates.MonthLocator()  # every month


gen = ((n.generators_t.p * n.generators_t.mu_upper)
       .groupby(n.generators.carrier, axis=1).sum())

sto = ((n.storage_units_t.p * n.storage_units_t.mu_upper)
       .groupby(n.storage_units.carrier, axis=1).sum())

capex = pd.concat([gen, sto], axis=1)

fig, ax = plt.subplots(figsize = (8, 3))

capex.div(1e6).plot(ax=ax)

ax.set_ylabel('Million €')
ax.set_xlabel('')
ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
fig.autofmt_xdate()


