#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 21:56:19 2020

@author: fabian
"""

import xarray as xr
import matplotlib.pyplot as plt
from config import color, to_total_symbol

plt.rc('font', family='serif')
plt.rc("text", usetex=False)

ca = xr.open_dataset(snakemake.input.costs)

fig, ax = plt.subplots()
ca.sum().to_array().to_series().reindex_like(color).dropna()\
    .rename(to_total_symbol).rename('')\
    .plot.pie(colors=color, explode=[.01]*len(ca), autopct='%1.0f%%')

fig.savefig(snakemake.output[0])