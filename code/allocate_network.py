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

from numpy.testing import assert_allclose
from pypsa.descriptors import nominal_attrs

from netallocation.utils import (reindex_by_bus_carrier as by_bus_carrier,
                                 get_as_dense_by_bus_carrier)
from netallocation.breakdown import (expand_by_sink_type,
                                     expand_by_source_type)
from netallocation.convert import peer_to_peer, virtual_patterns
from netallocation.cost import (nodal_co2_price, snapshot_weightings)
from config import source_dims
# from helpers import adjust_shadowprice

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('allocate_network', nname='de50')


n = pypsa.Network(snakemake.input.network)

aggregated = snakemake.config['aggregated_coupling']
method = snakemake.config['method']

if 'co2_limit' in n.global_constraints.index:
    co2_price = 0
else:
    co2_price = snakemake.config['costs']['emission_prices']['co2']
if 'lv_limit' in n.global_constraints.index:
    for c in n.branch_components:
        if n.df(c).empty: continue
        mu_upper = (n.global_constraints.at['lv_limit', 'mu'] * n.df(c).length)
        nom = nominal_attrs[c]
        n.df(c)['mu_upper_'+nom] += mu_upper


ds = ntl.allocate_flow(n, method=method, aggregated=aggregated)  # q=0

if snakemake.config['alloc_to_load_only']:
    ds = expand_by_sink_type(ds.chunk(dict(snapshot=5)), n)
    ds = ds.sel(sink_carrier='Load', drop=True)

A = expand_by_source_type(peer_to_peer(ds, n).peer_to_peer, n, chunksize=5)
A_f = virtual_patterns(ds, n, q=0).virtual_flow_pattern


# emission
ep = nodal_co2_price(n, price=co2_price) * snapshot_weightings(n)
ep = ep.rename(source_dims)
C_emission = A * ep


# opex
comps = ['Generator', 'StorageUnit', 'Store']
o = get_as_dense_by_bus_carrier(n, 'marginal_cost', comps).rename(source_dims)
o = o * snapshot_weightings(n) - ep.reindex_like(o, fill_value=0)
O = A * o


# =============================================================================
# WATCH OUT! In the paper the constraints always have a '<=' sign which lead
# to positive shadowprices.
# This is not the case for the pypsa.linopf implementation -> we have
# to adjust the mu's in some cases in order to take into account that they're
# not positive.
# =============================================================================

def scarcity_share(c):
    """Get the ratio of payments allocated to the scarcity constraint."""
    mu_upper = 'mu_upper_' + nominal_attrs[c]
    return (- n.df(c)[mu_upper] /
            (n.df(c).capital_cost - n.df(c)[mu_upper]))


# capex generator
c = 'Generator'
mu_gen = by_bus_carrier(n.pnl(c).mu_upper, c, n).rename(source_dims)
C_gen = mu_gen * A.sel(source_carrier=n.df(c).carrier.unique())
C_sca_gen = C_gen * by_bus_carrier(scarcity_share(c), c, n).rename(source_dims)
# rename for not conflicting later
C_gen = C_gen.rename(source_carrier='source_carrier_gen')

if not 'test' in snakemake.input.network:
    scarcity = n.df(c).p_nom_max.replace(np.inf, 0) @ n.df(c).mu_upper_p_nom
    subsidy = n.df(c).p_nom_min @ n.df(c).mu_lower_p_nom

    investment = n.df(c).p_nom_opt @ n.df(c).capital_cost
    revenue = (C_gen.sum() + subsidy + scarcity)
    assert_allclose(investment, revenue, rtol=1e-2, atol=np.inf)


# capex storage
c = 'StorageUnit'
if not n.df(c).empty:
    mu_sus = (n.pnl(c).mu_state_of_charge / n.df(c).efficiency_dispatch
              + n.pnl(c).mu_upper_p_dispatch + n.pnl(c).mu_lower_p_dispatch)
    mu_sus = by_bus_carrier(mu_sus, c, n).rename(source_dims)
    C_sus = mu_sus * A.sel(source_carrier=n.df(c).carrier.unique())
    C_sus = C_sus.rename(source_carrier='source_carrier_sus')



# capex branch
names = ['component', 'branch_i']
comps = set(A_f.component.data)
mu = pd.concat([n.pnl(c).mu_upper + n.pnl(c).mu_lower for c in comps],
               axis=1, names=names, keys=comps)
mu = xr.DataArray(mu, dims=['snapshot', 'branch'])
C_branch = (A_f * mu).rename(bus='sink')

share_scarcity = pd.concat([scarcity_share(c) for c in comps], keys=comps,
                           names=names)
share_scarcity = xr.DataArray(share_scarcity, dims='branch')
C_sca_branch = C_branch * share_scarcity

# collect cost allocations
ca = xr.Dataset({'one_port_operational_cost': O,
                 'co2_cost': C_emission,
                 'generator_investment_cost': C_gen,
                 'storage_investment_cost': C_sus,
                 'branch_investment_cost': C_branch})

sa = xr.Dataset({'generator_scarcity_cost': C_sca_gen,
                 'branch_scarcity_cost': C_sca_branch})

payments = ca.sum(['source', 'branch'])


ca.reset_index('branch').to_netcdf(snakemake.output.costs)
sa.reset_index('branch').to_netcdf(snakemake.output.scarcity)
payments.to_netcdf(snakemake.output.payments)
A.sum(['source']).mean('snapshot').to_netcdf(snakemake.output.power_mix)

# the total revenue for generators and branches is obtained by
# ca.sum() payed by loads and storage_units

# If would be possible to track the redistribution of payments by the storages.
# However, this would require knowledge about when the power supplied by a storage
# was initially charged.


