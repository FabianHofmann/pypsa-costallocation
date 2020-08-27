#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 12:53:08 2020

@author: fabian
"""

import pypsa
import geopandas as gpd
import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib as mpl
import pypsa
from config import color, to_symbol, sink_dims, to_explanation
import netallocation as ntl
from matplotlib.pyplot import Line2D
from pypsa.plot import projected_area_factor
from netallocation.plot_helpers import (make_handler_map_to_scale_circles_as_in,
                                        make_legend_circles_for)


if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('plot_allocated_payment', nname='de50bf',
                               method='ptpf', power='net', sink='DE0 1')


n = pypsa.Network(snakemake.input.network)
costs = xr.open_dataset(snakemake.input.costs)
sink = snakemake.wildcards.sink
regions = gpd.read_file(snakemake.input.regions).set_index('name')

if sink == 'lowest-lmp':
    sink = n.buses_t.marginal_price.mean().idxmin()
elif sink == 'highest-lmp':
    sink = n.buses_t.marginal_price.mean().idxmax()


payment = (xr.open_dataset(snakemake.input.costs)
           .sel(sink=sink, drop=True)
           .fillna(0)
           .mean('snapshot')
           .rename(source='bus', source_carrier='carrier')
           .set_index(branch=['component', 'branch_i']))


branch_widths = payment.branch_investment_cost.to_pandas()
branch_colors = (np.sign(payment.branch_investment_cost.to_pandas())
                 .replace({1: 'cadetblue', -1:'indianred'}))

emission_cost = payment[['co2_cost']].sum('carrier').to_dataframe().stack()
capex = pd.concat([payment.generator_investment_cost.to_series(),
                   payment.storage_investment_cost.to_series()]).sort_index()
opex = (payment.one_port_operational_cost.to_series()
        .rename(index=lambda s: s + ' opex', level=1))
bus_sizes = pd.concat([emission_cost, capex, opex]).sort_index()

capex_colors = n.carriers.color
opex_alpha = 1
opex_colors = (capex_colors
               .apply(lambda c: mpl.colors.to_rgba(c, opex_alpha))
                .apply(lambda c: mpl.colors.to_hex(c, True))
               .rename(lambda s: s + ' opex'))
emission_color = color.loc[['co2_cost']]
bus_colors = pd.concat([capex_colors, opex_colors, emission_color])


# %%

bus_scale = 1e-6
branch_scale = 3e-3

fig, ax = plt.subplots(figsize=(4.5,4.5), subplot_kw={'projection': ccrs.EqualEarth()})
n.plot(line_widths=branch_widths['Line'] * branch_scale,
       link_widths=branch_widths['Link'] * branch_scale,
       line_colors=branch_colors['Line'],
       link_colors=branch_colors['Link'],
       ax=ax,
       bus_alpha=None,
        geomap='10m',
       boundaries=regions.total_bounds[[0,2,1,3]],
       bus_sizes=bus_sizes * bus_scale,
       bus_colors=bus_colors,
       )
regions.loc[[sink]].plot(ax=ax, transform=ccrs.PlateCarree(), aspect='equal')


legend_colors = (n.carriers.set_index('nice_name').color
                 .append(emission_color.rename(to_explanation)))
fig.legend(
    *ntl.plot.handles_labels_for(legend_colors),
    loc="upper left",
    bbox_to_anchor=(1, 1),
    title='Carrier',
    frameon=False,
)



# legend generator capacities
reference_caps = [100e3, 10e3]
scale = 1 / bus_scale / projected_area_factor(ax)**2
handles = make_legend_circles_for(reference_caps, scale=scale,
                                  facecolor="w", edgecolor='grey',
                                  alpha=.5)
labels = ["%ik €"%(s / 1e3) for s in reference_caps]
handler_map = make_handler_map_to_scale_circles_as_in(ax)
legend = fig.legend(handles, labels,
                loc="upper left", bbox_to_anchor=(1., 0.4),
                frameon=False,
                handler_map=handler_map)
fig.add_artist(legend)

# legend transmission capacitues
handles, labels = [], []
reference_caps = [5,-2]
handles = [Line2D([0], [0], color=c, linewidth=abs(s)*1e3*branch_scale)
           for s, c in zip(reference_caps, ['cadetblue', 'indianred'])]
labels = ['%ik €'%s for s in reference_caps]

legend = fig.legend(handles, labels,
                    loc="lower left", bbox_to_anchor=(1, .1),
                    frameon=False,
                    )
fig.artists.append(legend)



fig.canvas.draw()
fig.tight_layout()
fig.savefig(snakemake.output[0], bbox_inches='tight')
