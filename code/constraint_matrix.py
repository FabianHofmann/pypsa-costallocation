#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 11:00:12 2020

@author: fabian
"""

import pypsa
from pypsa.linopt import set_int_index
import pandas as pd
import netallocation as ntl
from scipy.sparse.linalg import inv as spinv
from scipy.sparse import diags
from numpy.linalg import inv
from xarray import DataArray
import numpy as np
from numpy.random import rand
from numpy.testing import assert_allclose
from pandas import DataFrame, Series
from numpy import array
from pypsa.linopt import get_var, get_con
from pandas import concat

# %% create example

n0 = '1'
d0 = 60
d1 = 131
o0 = 52
o1 = 101
o2 = 156
c0 = 1060
c1 = 1050
c2 = 1020
# c_l = 100

n = pypsa.Network()

# old example
# n.madd('Bus', [n0, n1], x=[0, 1], y=[0, 0.4])
# n.madd('Load', [0, 1], bus=[n0, n1], p_set=[d0, d1])
# n.madd('Generator', [0, 1], bus=[n0, n1], p_nom_extendable=True,
#         marginal_cost=[o0, o1], capital_cost=[c0, c1], p_nom_max=100)
# n.madd('Line', ['1'], s_nom_extendable=True, x=0.01, bus0=[n0], bus1=[n1],
#         capital_cost=c_l)


# 1. example
# n.set_snapshots(list(range(2)))
# n.madd('Bus', [n0])
# n.madd('Load', [0], bus=[n0], p_set=pd.DataFrame({0: [d0, d1]}))
# n.madd('Generator', [0, 1, 2], bus=n0, p_nom_extendable=False,
#         marginal_cost=[o0, o1, o2], capital_cost=[c0, c1, c2], p_nom=[80, 51, 40],
#         p_min_pu=pd.DataFrame(rand(2,3)/100, range(2), range(3)),
#         p_nom_min = rand(3)/100,
#         carrier=list('abc'))

# 2. example
n.set_snapshots(list(range(2)))
n.madd('Bus', [n0])
n.madd('Load', [0], bus=[n0], p_set=pd.DataFrame({0: [d0, d1]}))
n.madd('Generator', [0, 1, 2], bus=n0, p_nom_extendable=True,
        marginal_cost=[o0, o1, o2], capital_cost=[c0, c1, c2],
        p_nom=[10, 11, 10],
        p_nom_max=100,
        p_min_pu=pd.DataFrame(rand(2,3)/100, range(2), range(3)),
        carrier=list('abc'))


# n = ntl.test.get_network_ac_dc()
# n.generators['p_nom'] = 0
# n.links['p_nom'] = 0
# n.lines['s_nom'] = 0
# n.set_snapshots(n.snapshots[[0]])


# %% add noise to all upper and lower bounds
for c, attr in pypsa.descriptors.nominal_attrs.items():
    if n.df(c).empty:
        continue

    n.df(c)[attr + '_min'] += rand(len(n.df(c)))/100
    n.df(c)[attr + '_max'] += rand(len(n.df(c)))/100

    min_pu = attr.replace('nom', 'min_pu')
    max_pu = attr.replace('nom', 'max_pu')

    default = n.df(c)[max_pu]
    n.pnl(c)[max_pu] = n.pnl(c)[max_pu]\
        .reindex(index=n.snapshots, columns=n.df(c).index).fillna(default)
    assert all(n.pnl(c)[max_pu].notnull())
    n.pnl(c)[max_pu] += rand(*n.pnl(c)[max_pu].shape)/100

    if min_pu in n.pnl(c):
        default = n.df(c)[min_pu]
        n.pnl(c)[min_pu] = n.pnl(c)[min_pu]\
            .reindex(index=n.snapshots, columns=n.df(c).index).fillna(default)
        assert all(n.pnl(c)[min_pu].notnull())
        n.pnl(c)[min_pu] += rand(*n.pnl(c)[min_pu].shape)/100


n.lopf(pyomo=False, keep_shadowprices=True, keep_references=True,
       solver_name='gurobi')


# %% get indeces of all variables and constraints
def to_series(df):
    """Standardize time and static variables."""
    return df.stack() if df.ndim > 1 else concat([df], keys=['static'])

vars, cons = {}, {}
for comp, attr in n.variables.index:
    vars[comp, attr] = to_series(get_var(n, comp, attr)).swaplevel(-1, -2)
vars = concat(vars, axis=0, names=['component', 'name', 'index', 'time'])
grb_varindex = Series({v.VarName: v.index for v in n.model.getVars()})\
                      .pipe(set_int_index)
vars_index = vars.map(grb_varindex).sort_values().index

for comp, attr in n.constraints.index:
    cons[comp, attr] = to_series(get_con(n, comp, attr)).swaplevel(-1, -2)
cons = concat(cons, axis=0, names=['component', 'name', 'index', 'time'])
grb_conindex = Series({c.ConstrName: c.index for c in n.model.getConstrs()})\
                      .pipe(set_int_index)
cons_index = cons.map(grb_conindex).sort_values().index

demand = (n.loads_t.p.groupby(n.loads.bus, axis=1)
          .sum().reindex(columns=n.buses.index, fill_value=0)
          .stack())


# %% extract allocations

tolerance = dict(atol=1e-5, rtol=1e-5)

O = n.model.getObjective()
A = n.model.getA().T # first dimension {vars} second {cons}
x = array([v.x for v in n.model.getVars()])
d = array([c.RHS  for c in n.model.getConstrs()])
c = np.zeros_like(x)
for t in range(O.size()): c[O.getVar(t).index] = O.getCoeff(t)
m = array([c.Pi for c in n.model.getConstrs()])

# remove the objective constant variable
nonzero_b = A.getnnz(1) > 0
A = A[nonzero_b]
x = x[nonzero_b]
c = c[nonzero_b]
B = (x * A).round(8) == d
N, M = A.shape

assert sum(B) == N
assert_allclose(x * A[:, B], d[B], **tolerance)

A_inv = spinv(A[:, B].tocsc()).T
lmp_i = Series(range(N), cons_index[B]).Bus.marginal_price
lmp = m[B][lmp_i]
r = A_inv[:, lmp_i] * diags(d[B][lmp_i])

assert_allclose(d[B][lmp_i], demand.values, **tolerance)
assert_allclose(x, A_inv * d[B], **tolerance)
assert_allclose(c * r, lmp * demand.values, rtol=1e-7, atol=1e-2)

r = DataFrame.sparse.from_spmatrix(r, vars_index, lmp_i.index).fillna(0)\
        .rename_axis(columns=['bus', 'snapshot'])


capex = r.loc['Generator', 'p_nom', :, 'static']
capex = capex.loc[:, (capex >= 0.1).any()]
opex = r.loc['Generator', 'p', :, :]


# # Now convert to indexed DataFrames
# A_ = DataFrame.sparse.from_spmatrix(A[:, B], vars_index, cons_index[B])
# A = DataFrame.sparse.from_spmatrix(A, vars_index, cons_index)
# A_inv = DataFrame.sparse.from_spmatrix(A_inv, vars_index, cons_index[B])
# D = DataFrame.from_spmatrix(D, vars_index, lmp_i.index)
# c = Series(c, vars_index)
# x = Series(x, vars_index)
# d = Series(d, cons_index)
# m = Series(m, cons_index)

# assert_allclose ((A_ @ A_inv).values, np.eye(N), **tolerance)
# assert_allclose(x @ A_, d[B], **tolerance)
# assert_allclose(x, d[B] @ A_inv, **tolerance)

# D = A_inv.reindex(['marginal_price'], level='name')
# demand = (n.loads_t.p.groupby(n.loads.bus, axis=1)
#           .sum().reindex(columns=n.buses.index, fill_value=0)
#           .stack().swaplevel(-1, -2))
# demand = concat({('Bus', 'marginal_price'): demand}, names=D.index.names)

# r = D.mul(demand, axis=0)
# lmp = m.reindex_like(demand)

# assert_allclose( r @ c / demand, lmp, **tolerance)

# lmp_decomposed = (r * c).T / demand

# xarray ansatz

# # convert them to xarray to keep track of the indeces
# coords = {'variables': vars.index, 'constraints': cons.index}
# A = DataArray(A.todense(), coords=coords, dims=coords)
# c = DataArray(c, coords={'variables': vars.index}, dims='variables')
# x = DataArray(x, coords={'variables': vars.index}, dims='variables')
# d = DataArray(d, coords={'constraints': cons.index}, dims='constraints')

# A_inv = A_.T.copy(data=inv(A_.data.to_scipy_sparse()))
# A_inv = A_.T.copy(data=np.linalg.inv(A_)) # coords have to be swapped
