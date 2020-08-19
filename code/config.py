#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:20:25 2020

@author: fabian
"""

import pandas as pd
import matplotlib.pyplot as plt

plt.rc('font', family='serif')
plt.rcParams['savefig.dpi'] = 200

source_dims = dict(bus='source', carrier='source_carrier')
sink_dims = dict(bus='sink', carrier='sink_carrier')
source_dims_r = dict(source='bus', source_carrier='carrier')
sink_dims_r = dict(sink='bus', sink_carrier='carrier')

source_carrier = ['source_carrier', 'source_carrier_gen', 'source_carrier_sus']


symbol = dict(generator_investment_cost = '\mathcal{C}^{G}',
              storage_investment_cost = '\mathcal{C}^{E}',
              branch_investment_cost = '\mathcal{C}^{F}',
              one_port_operational_cost = '\mathcal{O}',
              one_port_investment_cost = '\mathcal{C}',
              branch_sparcity_cost = '\mathcal{R}^{F sparcity}',
              one_port_sparcity_cost = '\mathcal{G sparcity}',
              co2_cost = '\mathcal{E}',
              lmp = '\lambda',
              demand = 'd')

to_symbol = {k: v + '_{n,t}' for k,v in symbol.items()}
to_static_symbol = {k: v + '_{n}' for k,v in symbol.items()}
to_total_symbol = {k: v for k,v in symbol.items()}

for d in [to_symbol, to_static_symbol, to_total_symbol]:
    d['nodal_demand_cost'] = d['lmp'] + '\,' + d['demand']
    for k in d:
        d[k] = '$'+ d[k] + '$'

to_explanation = {'one_port_operational_cost': 'OPEX',
                  'one_port_investment_cost': 'Production & Storage CAPEX',
                  'co2_cost': 'Emission Cost',
                  'generator_investment_cost': 'Production CAPEX',
                  'storage_investment_cost': 'Storage CAPEX',
                  'branch_investment_cost': 'Transmission CAPEX',
                  'branch_sparcity_cost': 'Sparcity Cost Transmission',
                  'one_port_sparcity_cost': 'Sparcity Cost Production',
                  'nodal_demand_cost': 'Nodal Payment'}

color = pd.Series({'one_port_operational_cost': 'darkkhaki',
                   'one_port_investment_cost': 'palevioletred',
                   'co2_cost': 'tomato',
                   'generator_investment_cost': 'palevioletred',
                   'storage_investment_cost': 'royalblue',
                   'branch_investment_cost': 'mediumaquamarine',
                   'nodal_demand_cost': 'cadetblue',
                   'remaining_cost': 'lightsteelblue'}).sort_index()

latex_names = dict(p_nom_opt = '$G$', marginal_cost = '$o$', capital_cost='$c$',
                   bus='$n$')
