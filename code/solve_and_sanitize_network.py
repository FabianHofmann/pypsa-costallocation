#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 22:16:58 2020

@author: fabian
"""

import pypsa
import numpy as np
import shutil


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('solve_and_sanitize_german_network', clusters=50)


solver_opts = snakemake.config['solving']['solver']
solver_name = solver_opts.pop('name')

n = pypsa.Network(snakemake.input.network)

n.generators.p_nom_max.update(n.generators.query('not p_nom_extendable').p_nom)
n.generators['p_nom_extendable'] = True

n.lopf(pyomo=False, solver_name=solver_name, solver_options=solver_opts,
        keep_shadowprices=True)

n.mremove('Line', n.lines.query('num_parallel == 0').index)

n.export_to_netcdf(snakemake.output.network)

shutil.copy(snakemake.input.regions, snakemake.output.regions)
