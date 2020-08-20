#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 13:26:32 2020

@author: fabian
"""

import pypsa
from pypsa.linopt import set_int_index, get_var, get_con, Dict
import pandas as pd
from scipy.sparse.linalg import inv as spinv
from scipy.sparse import diags
import numpy as np
from numpy.random import rand
from numpy.testing import assert_allclose
from pandas import DataFrame, Series
from numpy import array
from pypsa.descriptors import nominal_attrs
from pandas import concat
import xarray as xr
import matplotlib as mpl

fmt = mpl.ticker.ScalarFormatter(useMathText=True)
fmt.set_powerlimits((0, 0))

def combine_oneports(ds):
    gen = ds.generator_investment_cost.rename(source_carrier_gen='source_carrier')
    sus = ds.storage_investment_cost.rename(source_carrier_sus='source_carrier')
    return ds.drop(['generator_investment_cost', 'storage_investment_cost'])\
            .drop_dims(['source_carrier_gen', 'source_carrier_sus']) \
            .assign(one_port_investment_cost=xr.concat([gen, sus], dim='source_carrier'))

def adjust_shadowprice(gamma, c, n, s=None):
    nom = nominal_attrs[c]
    upper = (n.df(c).capital_cost /  (n.df(c).capital_cost - n.df(c)['mu_upper_'+nom]))
    if c == 'StorageUnit':
        p = 'p_dispatch'
    elif c in n.branch_components:
        p = 'p0'
    else:
        p = 'p'
    lower = (n.df(c)['mu_lower_'+nom] * n.df(c)[nom +'_opt'] / n.pnl(c)[p].sum())
    lower = lower.replace(np.inf, 0)

    return (gamma + lower) * upper

    # denom = (n.df(c).capital_cost - n.df(c)['mu_upper_'+nom] - n.df(c)['mu_lower_'+nom])
    # return gamma * (n.df(c).capital_cost /  denom)


