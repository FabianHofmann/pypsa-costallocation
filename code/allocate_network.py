#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 20:56:10 2020

@author: fabian
"""

import netallocation as ntl
import pypsa
from netallocation.utils import reindex_by_bus_carrier as by_bus_carrier

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('allocate_network', nname='de50',
                                   method='ptpf', power='net')

n = pypsa.Network(snakemake.input.network)
# n.set_snapshots(n.snapshots[:10])
method = snakemake.wildcards.method
aggregated = snakemake.wildcards.power == 'net'
co2_price = snakemake.config['costs']['emission_prices']['co2']


ds = ntl.allocate_flow(n, method=method, aggregated=aggregated, dask=False)  # q=0
ds = ds.chunk(dict(snapshot=5))
ca = ntl.allocate_cost(n, method=ds, chunksize=5, co2_price=co2_price)
ca['one_port_operational_cost'] = ca['one_port_operational_cost'] - \
                                  ca['co2_cost']

payments = ca.sum([d for d in ca.dims if d not in ["payer", "snapshot"]])

ca.reset_index('receiver_transmission_cost').to_netcdf(snakemake.output.costs)
payments.to_netcdf(snakemake.output.payments)


# pr = nodal_production_revenue(n).rename(bus="payer")
# dc = nodal_demand_cost(n).rename(bus="payer")
