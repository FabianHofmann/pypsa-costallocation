#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 10:39:11 2020

@author: fabian
"""

import pypsa
import numpy as np
import matplotlib.pyplot as plt
import netallocation as ntl
import cartopy.crs as ccrs
import geopandas as gpd


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('plot_german_network', clusters=50)


n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions)

fig, ax = plt.subplots(subplot_kw={"projection": ccrs.EqualEarth()},
                       figsize=(6, 5))
bus_sizes = n.generators.groupby(["bus", "carrier"]).p_nom_opt.sum()
plot = n.plot(
    bus_sizes=bus_sizes / 2e5,
    line_widths=n.lines.s_nom_opt / 3e3,
    link_widths=n.links.p_nom_opt / 3e3,
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
fig.savefig(snakemake.output[0])



