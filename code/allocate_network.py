#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 20:56:10 2020

@author: fabian
"""

import netallocation as ntl
import pypsa

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_german_network', clusters=50,
                                   method='ptpf', power='net')

n = pypsa.Network(snakemake.input.network)
# n.set_snapshots(n.snapshots[:10])
method = snakemake.wildcards.method
aggregated = snakemake.wildcards.power == 'net'

ds = ntl.allocate_flow(n, method=method, aggregated=aggregated,)  # q=0
ca = ntl.allocate_cost(n, method=ds)
payments = ca.sum([d for d in ca.dims if d not in ["payer", "snapshot"]])

ca.reset_index('receiver_transmission_cost').to_netcdf(snakemake.output.costs)
payments.to_netcdf(snakemake.output.payments)


# pr = nodal_production_revenue(n).rename(bus="payer")
# dc = nodal_demand_cost(n).rename(bus="payer")
