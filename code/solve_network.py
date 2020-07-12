#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 10:39:11 2020

@author: fabian
"""

import pypsa
import numpy as np

n = pypsa.Network(snakamake.input)

fig, ax = plt.subplots(subplot_kw={"projection": ccrs.PlateCarree()},
                       figsize=(6, 5))
bus_sizes = n.generators.groupby(["bus", "carrier"]).p_nom_opt.sum()
plot = n.plot(
    bus_sizes=bus_sizes / 1e5,
    line_widths=n.lines.s_nom_opt / 1e3,
    link_widths=n.links.p_nom_opt / 1e3,
    ax=ax,
)
ax.legend(
    *ntl.plot.handles_labels_for(n.carriers.color),
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    ncol=1,
)
fig.canvas.draw()
fig.tight_layout()
fig.savefig(snakemake.output)



# mu_co2 = 100


# n.calculate_dependent_values()
# n.generators.p_nom_max.update(n.generators.query('not p_nom_extendable').p_nom)
# n.generators["p_nom_extendable"] = True
# n.lines["s_nom_min"] = 0
# n.mremove("StorageUnit", n.storage_units.index)
# n.mremove("Line", n.lines.query("x_pu == inf").index)
# n.remove("Bus", "DE0 0")

# em = n.generators.carrier.map(n.carriers.co2_emissions) / \
#         n.generators.efficiency
# n.global_constraints = n.global_constraints.drop('CO2Limit')
# n.generators['marginal_cost'] += em * mu_co2


# n.lopf(pyomo=False, keep_shadowprices=True, solver_name="gurobi",
#        solver_options={'crossover': 0, 'method': 2})

# # reinsert mu
# n.global_constraints.at['CO2Limit', 'mu'] = mu_co2
# n.export_to_netcdf(snakamake.output)
