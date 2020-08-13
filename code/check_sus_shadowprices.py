#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 11:06:12 2020

@author: fabian
"""

import pypsa


n = pypsa.Network('/home/fabian/papers/costallocation/resources/de10gf.nc')

close = lambda df1, df2: ((df1/df2).fillna(1).round(2) == 1).all()

# Generators

g = n.generators
gt = n.generators_t

assert close(g.capital_cost, (gt.mu_upper.mul(gt.p_max_pu, fill_value=1)).sum()
                              + g.mu_upper_p_nom + g.mu_lower_p_nom)

g.eval('p_nom_opt * capital_cost')


# Lines/Links




# Storage Units

s = n.storage_units
st = n.storage_units_t


# This must be equal
assert close(s.capital_cost, (st.mu_upper_p_dispatch + st.mu_upper_p_store +
                              st.mu_upper_state_of_charge * s.max_hours).sum())


assert close((st.mu_lower_p_dispatch + st.mu_upper_p_dispatch + st.mu_state_of_charge / s.efficiency_dispatch)['DE0 0 battery'], n.buses_t.marginal_price['DE0 0'])


