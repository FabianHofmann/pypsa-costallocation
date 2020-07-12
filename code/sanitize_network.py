#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 22:16:58 2020

@author: fabian
"""

import pypsa
import numpy as np


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('sanitize_network', clusters=50)




n = pypsa.Network(snakemake.input.network)

n.mremove('Line', n.lines.query('num_parallel == 0').index)
n.add("GlobalConstraint", "CO2Limit_artificial",
      carrier_attribute="co2_emissions", sense="<=",
      constant=np.inf, mu=snakemake.config['costs']['emission_prices']['co2'])

n.export_to_netcdf(snakemake.output.network)

