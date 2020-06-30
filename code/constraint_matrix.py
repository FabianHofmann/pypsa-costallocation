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

n = networks.n2_t1_g2_w()
noisy_lopf(n)

print(n.generators_t.p)

s = get_linear_system(n)
A_, A_inv, c, x, d, m, r = (s[k].round(2) for k in
                            ['A_', 'A_inv', 'c', 'x', 'd', 'm', 'r'])

print(r)


C = r.mul(c, 0)

def pinv(df):
    return pd.DataFrame(np.linalg.pinv(df), *df.T.axes).round(3)

def inv(df):
    return pd.DataFrame(np.linalg.inv(df), *df.T.axes).round(3)