#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 23 17:33:26 2020

@author: fabian
"""

import matplotlib.pyplot as plt
import pandas as pd
from config import color
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.axisartist.axislines import SubplotZero

plt.rc('patch', facecolor='k')


cost_i = ['CAPEX Production', 'CAPEX Transmission',
          'OPEX Production', 'Other Cost Terms \n(CO$_2$, Maintenance etc.)']

lmp = pd.DataFrame([[25, 8, 15, 5]], index=['lmp'], columns=cost_i)
c = pd.Series(color[['one_port_investment_cost', 'branch_investment_cost',
                     'one_port_operational_cost', 'co2_cost']].values,
              index=cost_i)




fig = plt.figure(figsize=(1.5,4))
ax = SubplotZero(fig, 111)
fig.add_subplot(ax)
ax.axis['yzero'].set_axisline_style("-|>")
ax.axis['yzero'].set_visible(True)
for direction in ["left", "right", "bottom", "top"]:
    # hides borders
    ax.axis[direction].set_visible(False)

# fig, ax = plt.subplots(figsize=(1.5,4))
# for edge in ['bottom', 'right', 'top']:
#     ax.spines[edge].set_visible(False)

lmp.plot.bar(stacked=True, ax=ax, legend=False, color=c, position=-0.3)
ax.set_xlim(left=0)

bottom = 0
for i, c in enumerate(lmp):
    y = bottom + .5 * lmp.at['lmp', c]
    # ax.bar(1, lmp.at['lmp', c], 0.2, bottom)
    bottom += lmp.at['lmp', c]
    ax.annotate(c, xy=(.67, y), xytext=(.9, y), zorder=8,
                arrowprops=dict(arrowstyle='-'),
                verticalalignment='center')

ax.set_ylabel('$\lambda_{n,t}$  [€/MWh]')
ax.set_xticks([])
ax.set_yticklabels([])


fig.tight_layout()
fig.savefig('../figures/price_decomposition.png', bbox_inches='tight')
