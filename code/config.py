#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:20:25 2020

@author: fabian
"""

import pandas as pd

to_symbol = dict(one_port_investment_cost = '$\mathcal{C}^{G}_{n,t}$',
                      branch_investment_cost = '$\mathcal{C}^{F}_{n,t}$',
                      one_port_operational_cost = '$\mathcal{O}_{n,t}$',
                      co2_cost = '$\mathcal{E}_{n,t}$',
                      nodal_demand_cost = '$\lambda_{n,t} \, d_{n,t}$')

to_static_symbol = dict(one_port_investment_cost = '$\mathcal{C}^{G}_{n}$',
                      branch_investment_cost = '$\mathcal{C}^{F}_{n}$',
                      one_port_operational_cost = '$\mathcal{O}_{n}$',
                      co2_cost = '$\mathcal{E}_{n}$',
                      nodal_demand_cost = '$\lambda_{n} \, d_{n}$')


color = pd.Series({'one_port_operational_cost': 'darkkhaki',
                   'co2_cost': 'tomato',
                   'one_port_investment_cost': 'palevioletred',
                   'branch_investment_cost': 'mediumaquamarine',
                   'nodal_demand_cost': 'cadetblue'})

latex_names = dict(p_nom_opt = '$G$', marginal_cost = '$o$', capital_cost='$c$',
                   bus='$n$')