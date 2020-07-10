#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:00:12 2020

@author: fabian
"""

import pypsa
import numpy as np
import pandas as pd
from helpers import get_linear_system, noisy_lopf
import networks

n = networks.n_ac_dc()
em = n.generators.carrier.map(n.carriers.co2_emissions) / \
        n.generators.efficiency
noisy_lopf(n)
s = get_linear_system(n)
A_, A_inv, c, x, d, m, r, C = (s[k].round(5) for k in
                            ['A_', 'A_inv', 'c', 'x', 'd', 'm', 'r', 'C'])
# C = C.loc[:, d.round(1)!=0]


sn = n.snapshots[0]
bus = 'London'

gen = n.generators_t.p.T.copy()
TC = n.objective

print(n.generators_t.p @ em)
print(m)


mu_co2 = n.global_constraints.mu.co2_limit
n.global_constraints = n.global_constraints.drop('co2_limit')

n.generators['marginal_cost'] += em * mu_co2

n.lopf(pyomo=False, keep_shadowprices=True, keep_references=True,
        solver_name='gurobi')
s = get_linear_system(n)
A_, A_inv, c, x, d, m, r, C = (s[k].round(5) for k in
                            ['A_', 'A_inv', 'c', 'x', 'd', 'm', 'r', 'C'])


# print(n.generators_t.p @ em)


# n.loads_t.p_set.loc[sn,bus] += 1

# n.lopf(pyomo=False, keep_shadowprices=True, keep_references=True,
#         solver_name='gurobi')
# print(n.generators_t.p.T - gen)
# print(n.objective - TC)


# n.global_constraints = n.global_constraints.drop('co2_limit')

# print('Objective diff:', objective - n.objective)
# print('CO2 cost diff:', total_co2_cost)