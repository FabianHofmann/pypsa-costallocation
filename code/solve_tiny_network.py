#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 17:50:06 2020

@author: fabian
"""

import pypsa
import netallocation as ntl

if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('solve_tiny_network')


n = ntl.test.get_network_ac_dc()
n.generators['p_nom_min'] = 500
n.generators['p_nom_max'] = 3000

# n.carriers = n.carriers.drop('battery')
n.add('StorageUnit', 'battery', bus='Frankfurt', capital_cost=50,
      p_nom_extendable=True, marginal_cost=10, carrier='battery')

n.lopf(pyomo=False, keep_shadowprices=True, solver_name='gurobi')

n.export_to_netcdf(snakemake.output[0])