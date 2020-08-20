#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 15:02:28 2020

@author: fabian
"""

import networks

n = networks.n_ac_dc()

n.lopf()
print(n.generators_t.p)
print(n.objective)

n.loads_t.p_set.iloc[0,0] += 1

n.lopf()
print(n.generators_t.p)
print(n.objective)
