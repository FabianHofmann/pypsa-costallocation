#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 12:26:49 2021

@author: fabian
"""

import netallocation as ntl
import pypsa

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('allocate_network', nname='test-de10bf',
                                   method='ptpf', power='net')


n = pypsa.Network(snakemake.input.network)


ds = ntl.allocate_revenue(n)
payments = ds.sum(['source', 'branch'])


ds.reset_index('branch').to_netcdf(snakemake.output.revenue)
payments.to_netcdf(snakemake.output.payments)



