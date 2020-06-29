#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 13:26:32 2020

@author: fabian
"""

import pypsa
from pypsa.linopt import set_int_index
import pandas as pd
from scipy.sparse.linalg import inv as spinv
from scipy.sparse import diags
import numpy as np
from numpy.random import rand
from numpy.testing import assert_allclose
from pandas import DataFrame, Series
from numpy import array
from pypsa.linopt import get_var, get_con, Dict
from pandas import concat



def add_noise(n, noise_order=1e-4):
    for c, attr in pypsa.descriptors.nominal_attrs.items():
        if n.df(c).empty:
            continue

        n.df(c)[attr + '_min'] += rand(len(n.df(c))) * noise_order
        n.df(c)[attr + '_max'] += rand(len(n.df(c))) * noise_order

        min_pu = attr.replace('nom', 'min_pu')
        max_pu = attr.replace('nom', 'max_pu')

        default = n.df(c)[max_pu]
        n.pnl(c)[max_pu] = n.pnl(c)[max_pu]\
            .reindex(index=n.snapshots, columns=n.df(c).index).fillna(default)
        assert all(n.pnl(c)[max_pu].notnull())
        n.pnl(c)[max_pu] += rand(*n.pnl(c)[max_pu].shape) * noise_order

        if min_pu in n.pnl(c):
            default = n.df(c)[min_pu]
            n.pnl(c)[min_pu] = n.pnl(c)[min_pu]\
                .reindex(index=n.snapshots, columns=n.df(c).index).fillna(default)
            assert all(n.pnl(c)[min_pu].notnull())
            n.pnl(c)[min_pu] += rand(*n.pnl(c)[min_pu].shape) * noise_order


def noisy_lopf(n, noise_order=1e-4):
    add_noise(n, noise_order)
    n.lopf(pyomo=False, keep_shadowprices=True, keep_references=True,
       solver_name='gurobi')

def to_series(df):
    """Standardize time and static variables."""
    return df.stack() if df.ndim > 1 else concat([df], keys=['static'])


def variable_index(n):
    vars = {}
    for comp, attr in n.variables.index:
        vars[comp, attr] = to_series(get_var(n, comp, attr)).swaplevel(-1, -2)
    vars = concat(vars, axis=0, names=['component', 'name', 'component_i', 'snapshot'])
    grb_varindex = Series({v.VarName: v.index for v in n.model.getVars()})\
                          .pipe(set_int_index)
    return vars.map(grb_varindex).sort_values().index

def constraint_index(n):
    cons = {}
    for comp, attr in n.constraints.index:
        cons[comp, attr] = to_series(get_con(n, comp, attr)).swaplevel(-1, -2)
    cons = concat(cons, axis=0, names=['component', 'name', 'component_i', 'snapshot'])
    grb_conindex = Series({c.ConstrName: c.index for c in n.model.getConstrs()})\
                          .pipe(set_int_index)
    return cons.map(grb_conindex).sort_values().index


def demand_per_bus(n):
    return (n.loads_t.p.groupby(n.loads.bus, axis=1).sum()
            .reindex(columns=n.buses.index, fill_value=0).stack())



def get_linear_system(n):
    tolerance = dict(atol=1e-5, rtol=1e-5)
    model = n.model

    O = model.getObjective()
    A = model.getA().T # first dimension {vars} second {cons}
    x = array([v.x for v in model.getVars()])
    d = array([c.RHS  for c in model.getConstrs()])
    c = np.zeros_like(x)
    for t in range(O.size()): c[O.getVar(t).index] = O.getCoeff(t)
    m = array([c.Pi for c in model.getConstrs()])

    # remove the objective constant variable
    nonzero_b = A.getnnz(1) > 0
    A = A[nonzero_b]
    x = x[nonzero_b]
    c = c[nonzero_b]
    B = (x * A).round(8) == d
    N, M = A.shape
    A_ = A[:, B]

    assert sum(B) == N
    assert_allclose(x * A_, d[B], **tolerance)

    A_inv = spinv(A_.tocsc()).T

    vars_index = variable_index(n)
    cons_index = constraint_index(n)
    s = Dict()

    s.A_ = DataFrame.sparse.from_spmatrix(A_, vars_index, cons_index[B])
    # A = DataFrame.sparse.from_spmatrix(A, vars_index, cons_index)
    s.A_inv = DataFrame.sparse.from_spmatrix(A_inv, vars_index, cons_index[B])
    s.c = Series(c, vars_index)
    s.x = Series(x, vars_index)
    s.d = Series(d[B], cons_index[B])
    s.m = Series(m[B], cons_index[B])

    cost_effective_i = s.d[lambda ds: ds!=0].index
    s.r = s.A_inv.loc[:, cost_effective_i] * s.d[cost_effective_i]

    assert_allclose(s.x, s.r.sum(1), **tolerance)

    return s

