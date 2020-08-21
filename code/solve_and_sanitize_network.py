#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 22:16:58 2020

@author: fabian
"""

import pypsa
import pandas as pd
import numpy as np
import shutil


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('solve_and_sanitize_german_network',
                                   clusters=10, field='bf')


solver_opts = snakemake.config['solving']['solver']
solver_name = solver_opts.pop('name')

n = pypsa.Network(snakemake.input.network)

n.mremove('StorageUnit', n.storage_units.query('carrier=="hydro"').index)
n.remove('Carrier', 'hydro')
for c in ['Generator', 'StorageUnit']:
    n.df(c).p_nom_max.update(n.df(c).query('not p_nom_extendable').p_nom)
    n.df(c)['p_nom_extendable'] = True

n.mremove('Line', n.lines.query('num_parallel == 0').index)

# remove preexisting infrastructure from objective function
for c, attr in pypsa.descriptors.nominal_attrs.items():
    n.df(c)[attr] = 0

# set lv limit here, as long as #175 is not solved
if 'line_volume_limit' in snakemake.config:
    links_dc_b = n.links.carrier == 'DC' if not n.links.empty else pd.Series()
    lines_s_nom = n.lines.s_nom.where(n.lines.type == '',
                                      np.sqrt(3) * n.lines.num_parallel *
                                      n.lines.type.map(n.line_types.i_nom) *
                                      n.lines.bus0.map(n.buses.v_nom))
    total_line_volume = ((lines_s_nom * n.lines['length']).sum() +
                         n.links.loc[links_dc_b].eval('p_nom * length').sum())
    line_volume = snakemake.config['line_volume_limit'] * total_line_volume
    n.add('GlobalConstraint', 'lv_limit',
          type='transmission_volume_expansion_limit',
          sense='<=', constant=line_volume, carrier_attribute='AC, DC')


if snakemake.wildcards.field == 'bf':
    for c, attr in pypsa.descriptors.nominal_attrs.items():
        n.df(c)[attr + '_min'] = n.df(c)[attr]
else:
    for c, attr in pypsa.descriptors.nominal_attrs.items():
        n.df(c)[attr + '_min'] = 0
        n.df(c)[attr] = 0

n.lopf(pyomo=False, solver_name=solver_name, solver_options=solver_opts,
        keep_shadowprices=True)

n.name = snakemake.output.network.split('/')[-1].replace('.nc', '')

n.export_to_netcdf(snakemake.output.network)

shutil.copy(snakemake.input.regions, snakemake.output.regions)
