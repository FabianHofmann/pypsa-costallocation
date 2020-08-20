#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 12:37:18 2020

@author: fabian
"""

import gurobi
from numpy import array
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


model = gurobi.read('tinylp.lp')
model.optimize()

O = model.getObjective()
A = model.getA().T.todense() # first dimension {vars} second {cons}
x = array([v.x for v in model.getVars()])
d = array([c.RHS  for c in model.getConstrs()])
c = np.zeros_like(x)
for t in range(O.size()): c[O.getVar(t).index] = O.getCoeff(t)
m = array([c.Pi for c in model.getConstrs()])


model2 = gurobi.read('tinylp2.lp')
model2.optimize()

O2 = model2.getObjective()
A2 = model2.getA().T.todense() # first dimension {vars} second {cons}
x2 = array([v.x for v in model2.getVars()])
d2 = array([c.RHS  for c in model2.getConstrs()])
c2 = np.zeros_like(x)
for t in range(O2.size()): c2[O2.getVar(t).index] = O2.getCoeff(t)
m2 = array([c.Pi for c in model2.getConstrs()])

print(x, x2)


#%% plot
N = 10
fig, (ax1, ax2) = plt.subplots(1,2, sharey=True, figsize=(8,5))
# ax1.axis('off')
ax1.set_aspect('equal', adjustable='box')
ax1.set_xlim(0,N-1)
# ax2.axis('off')
ax2.set_xlim(0,N-1)
ax2.set_aspect('equal', adjustable='box')

f = r'$2 x_0+4 x_1$'
con1 = r'$ x_1 + x_2 = 7 $'
con2 = r'$x_1 + 0.5 x_2 \le 4$'
mu2 = r'$\perp \mu = -4 $'
fmod = r'$5 x_0 + 5  x_1$'


carray = pd.DataFrame([[c @ np.array([x1,x2]) for x1 in range(N)]
                       for x2 in range(N)])
xx = np.linspace(0,N,200)
y = pd.DataFrame(d - (A[0].T * xx).T, index=xx) / np.asarray(A[1])
binding_con = y[0]
upper_bound = y[3]

ax1.contourf(carray, levels=50)
binding_con[binding_con>=.0].plot(ax=ax1, color='grey', lw=2.5, ls='dashed')
ax1.fill_between(xx[upper_bound>=0], upper_bound[upper_bound>=0], alpha=0.2,
                color='red')
ax1.scatter(*x, marker='8', color='white', zorder=8)
ax1.text(3.7,3.7, con1, color='white')
ax1.text(.7,1.5, con2, color='white')
ax1.text(.7,1., mu2, color='white')
ax1.text(*(x + (.2, .2)), '$x^* = (1, 6)$', color='white')
ax1.set_title(r'$f = \sum_i c_i x_i = 2 x_0 + 4 x_1$')

carray2 = pd.DataFrame([[c2 @ np.array([x1,x2]) for x1 in range(N)]
                       for x2 in range(N)])
y = pd.DataFrame(d2 - (A2[0].T * xx).T, index=xx) / np.asarray(A2[1])
binding_con = y[0]

binding_con[binding_con>=.0].plot(ax=ax2, color='grey', lw=2.5, ls='dashed')
ax2.contourf(carray2, levels=50, vmin=carray.min().min(), vmax=carray.max().max())
ax2.scatter(*x2, marker='8', color='white', zorder=8)
ax2.set_xlim(0)

ax2.text(3.7,3.7, con1, color='white')
ax2.text(*(x2 + (.2, .2)), '$x^* = (0, 7)$', color='white')
ax2.set_title(r'$\tilde{f} = \sum_i \tilde{c}_i x_i ='
              ' (5 + \epsilon) x_0 + (5 + \epsilon / 2)  x_1$')

fig.tight_layout()
fig.savefig('../figures/shifted_objective.png')