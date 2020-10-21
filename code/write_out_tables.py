#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 13:21:30 2020

@author: fabian
"""
import pypsa
import xarray as xr
import pandas as pd
import os
from helpers import combine_oneports, load

if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('write_tables', nname='test-de10bf',
                               method='ptpf', power='net', sink='DE0 1')


n = pypsa.Network(snakemake.input.network)
costs = xr.open_dataset(snakemake.input.costs)

outdir = snakemake.output[0]
os.makedirs(outdir, exist_ok=True)

# %%
# capital and operational prices
cols = ['marginal_cost', 'capital_cost']
n.lines['marginal_cost'] = 0
adjust_cp = (lambda df: df.assign(capital_cost=df.capital_cost/df.length)
             if 'length' in df else df)
prices = pd.concat({c: adjust_cp(n.df(c)).groupby('carrier')[cols].mean()
           for c in ['Generator', 'StorageUnit', 'Line', 'Link']}, )

index = {'Link': 'Line', 'StorageUnit': 'Storage'}
columns = {'marginal_cost': r'o [\euro/MWh]', 'capital_cost': r'c [k\,\euro/MW]$^*$'}
prices['capital_cost'] /= 1e3
prices = (prices.rename(columns=columns).rename(index, level=0)
          .rename(n.carriers.nice_name, level=1)).round(3)

with open(snakemake.output[0] + '/prices.tex', 'w') as file:
    file.write(prices.replace(0., '').to_latex(escape=False))

# %%

# ds = costs.sum(['source', 'sink', 'snapshot', 'branch']) / load(n).sum(['snapshot', 'sink'])

