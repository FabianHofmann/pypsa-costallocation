#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 20:56:10 2020

@author: fabian
"""

import netallocation as ntl
import pypsa
import xarray as xr
import pandas as pd

from netallocation.utils import (reindex_by_bus_carrier as by_bus_carrier,
                                 get_as_dense_by_bus_carrier)
from netallocation.breakdown import (by_carriers as breakdown_carriers,
                                     expand_by_sink_type,
                                     expand_by_source_type)
from netallocation.convert import peer_to_peer, virtual_patterns
from netallocation.cost import (nodal_co2_price, snapshot_weightings,
                                allocate_one_port_operational_cost,
                                allocate_branch_investment_cost)
from config import source_dims, sink_dims

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('allocate_network', nname='acdc',
                                   method='ptpf', power='net')

n = pypsa.Network(snakemake.input.network)

method = snakemake.wildcards.method
aggregated = snakemake.wildcards.power == 'net'
if 'co2_limit' in n.global_constraints.index:
    co2_price = 0
else:
    co2_price = snakemake.config['costs']['emission_prices']['co2']


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


# capex one port
c = 'Generator'
mu_gen = by_bus_carrier(n.pnl(c).mu_upper, c, n)

c = 'StorageUnit'
if not n.df(c).empty:
    mu_sus = (n.pnl(c).mu_state_of_charge / n.df(c).efficiency_dispatch
              - n.pnl(c).mu_upper_p_dispatch)
    # mu_sus_ch = (n.pnl(c).mu_state_of_charge * n.df(c).efficiency_store
    #              - n.pnl(c).mu_lower_p_store)
    mu_sus = by_bus_carrier(mu_sus, c, n)
    mu = xr.concat([mu_gen, mu_sus], dim='carrier').rename(source_dims).fillna(0)
else:
    mu = mu_gen.rename(source_dims).fillna(0)
A_capex = mu * A


# capex branch
names = ['component', 'branch_i']
comps = set(A_f.component.data)
mu = pd.concat({c: n.pnl(c).mu_upper + n.pnl(c).mu_lower for c in comps},
               axis=1, names=names)
mu = xr.DataArray(mu, dims=['snapshot', 'branch'])
A_capex_branch = (A_f * mu).rename(bus='sink')


allocate_branch_investment_cost(A_f, n).rename(bus='sink')

ca = xr.Dataset({'one_port_operational_cost': A_opex,
                 'co2_cost': A_emission,
                 'one_port_investment_cost': A_capex,
                 'branch_investment_cost': A_capex_branch})
ca = expand_by_sink_type(ca, n)

payments = ca.sum(['source', 'source_carrier', 'branch'])
nodal_payments = payments.to_array().sum(['variable', 'sink_carrier'])


ca.reset_index('branch').to_netcdf(snakemake.output.costs)
payments.to_netcdf(snakemake.output.payments)



# pr = nodal_production_revenue(n).rename(bus="payer")
# dc = nodal_demand_cost(n).rename(bus="payer")
