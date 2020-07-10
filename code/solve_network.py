#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 10:39:11 2020

@author: fabian
"""

import pypsa
import numpy as np

name = 'elec_s_50_ec_lvopt_Co2L-3H'
n = pypsa.Network(f"../{name}.nc")
n.calculate_dependent_values()
n.generators["p_nom_extendable"] = True
n.lines["s_nom_min"] = 0
# n.lines['s_nom_max'] = 6000
# n.generators['p_nom_max'] = np.inf
n.mremove("StorageUnit", n.storage_units.index)
n.mremove("Line", n.lines.query("x_pu == inf").index)
n.remove("Bus", "DE0 0")
n.lopf(pyomo=False, keep_shadowprices=True, solver_name="gurobi",
       solver_options={'crossover': 0, 'method': 2})
n.to_netcdf('solved_germany_50.nc')