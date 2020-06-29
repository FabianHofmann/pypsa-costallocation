#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:00:12 2020

@author: fabian
"""

import pypsa
import pandas as pd
from helpers import get_linear_system, noisy_lopf
import networks

n = networks.n2_t1_g2_w()

noisy_lopf(n)
s = get_linear_system(n)
A_, A_inv, c, x, d, m, r = s.A_, s.A_inv, s.c, s.x, s.d, s.m, s.r
