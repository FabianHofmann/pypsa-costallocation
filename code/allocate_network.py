#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 20:56:10 2020

@author: fabian
"""

import netallocation as ntl
import pypsa
import xarray as xr

from netallocation.utils import (reindex_by_bus_carrier as by_bus_carrier,
                                 get_as_dense_by_bus_carrier)
from netallocation.breakdown import (by_carriers as breakdown_carriers,
                                     expand_by_sink_type,
                                     expand_by_source_type)
from netallocation.convert import peer_to_peer, virtual_patterns
from netallocation.cost import (nodal_co2_price, snapshot_weightings,
                                allocate_one_port_operational_cost,
                                allocate_branch_investment_cost)
from config import source_dims

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('allocate_network', nname='de10gf',
                                   method='ptpf', power='net')

n = pypsa.Network(snakemake.input.network)

method = snakemake.wildcards.method
aggregated = snakemake.wildcards.power == 'net'
co2_price = snakemake.config['costs']['emission_prices']['co2']


ds = ntl.allocate_flow(n, method=method, aggregated=aggregated, dask=False)  # q=0
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
o = o - ep.reindex_like(o, fill_value=0)
A_opex = A * o


# capex one port
c = 'Generator'
mu_gen = by_bus_carrier(n.pnl(c).mu_upper, c, n)

c = 'StorageUnit'
mu_sus = (- n.pnl(c).mu_state_of_charge / n.df(c).efficiency_dispatch -
          n.pnl(c).mu_lower_p_dispatch)
mu_sus = by_bus_carrier(mu_sus, c, n)
mu = xr.concat([mu_gen, mu_sus], dim='carrier').rename(source_dims).fillna(0)
A_capex = mu * A


# capex branch
A_capex_branch = allocate_branch_investment_cost(A_f, n).rename(bus='sink')

ca = xr.Dataset({'one_port_operational_cost': A_opex,
                 'co2_cost': A_emission,
                 'one_port_investment_cost': A_capex,
                 'branch_investment_cost': A_capex_branch})
ca = expand_by_sink_type(ca, n)

payments = ca.sum([d for d in ca.dims if d not in
                   ["sink", "sink_carrier", "snapshot"]])\
             .sel(sink_carrier='Load', drop=True)

ca.reset_index('branch').to_netcdf(snakemake.output.costs)
payments.to_netcdf(snakemake.output.payments)


# pr = nodal_production_revenue(n).rename(bus="payer")
# dc = nodal_demand_cost(n).rename(bus="payer")
