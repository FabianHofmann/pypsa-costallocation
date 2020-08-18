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

to_symbol = dict(one_port_investment_cost = '$\mathcal{C}^{G}_{n,t}$',
                 branch_investment_cost = '$\mathcal{C}^{F}_{n,t}$',
                 one_port_operational_cost = '$\mathcal{O}_{n,t}$',
                 branch_sparcity_cost = '$\mathcal{R}^{F sparcity}_{n,t}$',
                 one_port_sparcity_cost = '$\mathcal{G sparcity}_{n,t}$',
                 co2_cost = '$\mathcal{E}_{n,t}$',
                 nodal_demand_cost = '$\lambda_{n,t} \, d_{n,t}$'
                 )

to_static_symbol = dict(one_port_investment_cost = '$\mathcal{C}^{G}_{n}$',
                        branch_investment_cost = '$\mathcal{C}^{F}_{n}$',
                        one_port_operational_cost = '$\mathcal{O}_{n}$',
                        branch_sparcity_cost = '$\mathcal{R}^{F sparcity}_{n}$',
                        one_port_sparcity_cost = '$\mathcal{G sparcity}_{n}$',
                        co2_cost = '$\mathcal{E}_{n}$',
                        nodal_demand_cost = '$\lambda_{n} \, d_{n}$')


to_total_symbol = dict(one_port_investment_cost = '$\mathcal{C}^{G}$',
                      branch_investment_cost = '$\mathcal{C}^{F}$',
                      one_port_operational_cost = '$\mathcal{O}$',
                      branch_sparcity_cost = '$\mathcal{R}^{F sparcity}$',
                      one_port_sparcity_cost = '$\mathcal{G sparcity}$',
                      co2_cost = '$\mathcal{E}$',
                      nodal_demand_cost = '$\mathcal{TC}$',
                      remaining_cost = '$\mathcal{R}$')

to_explanation = {'one_port_operational_cost': 'OPEX',
                  'co2_cost': 'Emission Cost',
                  'one_port_investment_cost': 'Production CAPEX',
                  'branch_investment_cost': 'Transmission CAPEX',
                  'branch_sparcity_cost': 'Sparcity Cost Production',
                  'one_port_sparcity_cost': 'Sparcity Cost Transmission',
}

color = pd.Series({'one_port_operational_cost': 'darkkhaki',
                   'co2_cost': 'tomato',
                   'one_port_investment_cost': 'palevioletred',
                   'branch_investment_cost': 'mediumaquamarine',
                   'nodal_demand_cost': 'cadetblue',
                   'remaining_cost': 'lightsteelblue'}).sort_index()

latex_names = dict(p_nom_opt = '$G$', marginal_cost = '$o$', capital_cost='$c$',
                   bus='$n$')
