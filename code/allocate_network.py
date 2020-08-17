#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 20:56:10 2020

@author: fabian
"""

import netallocation as ntl
import pypsa
import xarray as xr
import numpy as np
import pandas as pd

from netallocation.utils import (reindex_by_bus_carrier as by_bus_carrier,
                                 get_as_dense_by_bus_carrier)
from netallocation.breakdown import (expand_by_sink_type,
                                     expand_by_source_type)
from netallocation.convert import peer_to_peer, virtual_patterns
from netallocation.cost import (nodal_co2_price, snapshot_weightings)
from config import source_dims
from helpers import adjust_shadowprice

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('allocate_network', nname='acdc',
                                   method='ptpf', power='net')


n = pypsa.Network(snakemake.input.network)

method = snakemake.wildcards.method.split('-')[0]
adjust_mu = '-adj' in snakemake.wildcards.method
aggregated = snakemake.wildcards.power == 'net'
if 'co2_limit' in n.global_constraints.index:
    co2_price = 0
else:
    co2_price = snakemake.config['costs']['emission_prices']['co2']
if 'lv_limit' in n.global_constraints.index:
    for c in n.branch_components:
        if n.df(c).empty: continue
        mu_upper = (n.global_constraints.at['lv_limit', 'mu'] * n.df(c).length)
        nom = pypsa.descriptors.nominal_attrs[c]
        n.df(c)['mu_upper_'+nom] += mu_upper


ds = ntl.allocate_flow(n, method=method, aggregated=aggregated)  # q=0
ds = ds.chunk(dict(snapshot=5))

A = expand_by_source_type(peer_to_peer(ds, n).peer_to_peer, n, chunksize=5)
A_f = virtual_patterns(ds, n, q=0).virtual_flow_pattern


# emission
ep = nodal_co2_price(n, price=co2_price) * snapshot_weightings(n)
ep = ep.rename(source_dims)
A_emission = A * ep


# opex
comps = ['Generator', 'StorageUnit', 'Store']
o = get_as_dense_by_bus_carrier(n, 'marginal_cost', comps).rename(source_dims)
o = o * snapshot_weightings(n) - ep.reindex_like(o, fill_value=0)
A_opex = A * o


# =============================================================================
# WATCH OUT! In the paper the constraints always have a '<=' sign which lead
# to positive shadowprices.
# This is not the case for the pypsa.linopf implementation -> we have
# to adjust the mu's in some cases in order to take into account that they're
# not positive.
# =============================================================================


# capex one port
c = 'Generator'
# calculate corrections for capacity restrictions
mu_gen = n.pnl(c).mu_upper
mu_gen = adjust_shadowprice(mu_gen, c, n) if adjust_mu else mu_gen
mu_gen = by_bus_carrier(mu_gen, c, n)


c = 'StorageUnit'
if not n.df(c).empty:
    # calculate corrections for capacity restrictions
    mu_sus = (n.pnl(c).mu_state_of_charge / n.df(c).efficiency_dispatch
              + n.pnl(c).mu_upper_p_dispatch + n.pnl(c).mu_lower_p_dispatch)
    mu_sus = adjust_shadowprice(mu_sus, c, n) if adjust_mu else mu_sus
    mu_sus = by_bus_carrier(mu_sus, c, n)
    mu = xr.concat([mu_gen, mu_sus], dim='carrier').rename(source_dims).fillna(0)
else:
    mu = mu_gen.rename(source_dims).fillna(0)
A_capex = mu * A


# capex branch
names = ['component', 'branch_i']
comps = set(A_f.component.data)
mu = []
for c in comps:
    mu_b = n.pnl(c).mu_upper + n.pnl(c).mu_lower
    mu.append(adjust_shadowprice(mu_b, c, n) if adjust_mu else mu_b)
mu = pd.concat(mu, axis=1, names=names, keys=comps)
mu = xr.DataArray(mu, dims=['snapshot', 'branch'])
A_capex_branch = (A_f * mu).rename(bus='sink')


# collect cost allocations
ca = xr.Dataset({'one_port_operational_cost': A_opex,
                 'co2_cost': A_emission,
                 'one_port_investment_cost': A_capex,
                 'branch_investment_cost': A_capex_branch})
ca = expand_by_sink_type(ca, n, chunksize=5)

payments = ca.sum(['source', 'source_carrier', 'branch'])
nodal_payments = payments.to_array().sum(['variable', 'sink_carrier'])


ca.reset_index('branch').to_netcdf(snakemake.output.costs)
payments.to_netcdf(snakemake.output.payments)

# the total revenue for generators and branches is obtained by
# ca.sum() payed by loads and storage_units

# If would be possible to track the redistribution of payments by the storages.
# However, this would require knowledge about when the power supplied by a storage
# was initially charged.


