#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 11:06:12 2020

@author: fabian
"""

import pypsa
from helpers import adjust_shadowprice

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('check_shadowprices', nname='de10bf')


n = pypsa.Network(snakemake.input[0])

tol = (n.objective + n.objective_constant)*1e-6

if 'lv_limit' in n.global_constraints.index:
    for c in n.branch_components:
        if n.df(c).empty: continue
        mu_upper = (n.global_constraints.at['lv_limit', 'mu'] * n.df(c).length)
        nom = pypsa.descriptors.nominal_attrs[c]
        n.df(c)['mu_upper_'+nom] += mu_upper

close = lambda df1, df2: ((df1 - df2).abs() < tol ).all()

# Generators

g = n.generators
gt = n.generators_t

assert close(g.capital_cost, (gt.mu_upper.mul(gt.p_max_pu, fill_value=1)).sum()
                              + g.mu_upper_p_nom + g.mu_lower_p_nom)

investment = g.eval('p_nom_opt * capital_cost')
revenue = (gt.mu_upper * gt.p).sum()

gamma = adjust_shadowprice(gt.mu_upper, 'Generator', n)
revenue = (gamma * gt.p).sum()
assert close(investment.where(gt.p.sum()>0, 0), revenue)



# Lines/Links

l = n.lines
lt = n.lines_t

investment = l.eval('s_nom_opt * capital_cost')
revenue = ((lt.mu_upper + lt.mu_lower) * lt.p0).sum()


gamma = lt.mu_upper + lt.mu_lower
gamma_adj = adjust_shadowprice(gamma, 'Line', n)
revenue = (gamma_adj * lt.p0).sum()
assert close(investment, revenue)


l = n.links
lt = n.links_t

investment = l.eval('p_nom_opt * capital_cost')
revenue = ((lt.mu_upper + lt.mu_lower) * lt.p0).sum()

lv_mu = - n.global_constraints.at['lv_limit', 'mu']
gamma = lt.mu_upper + lt.mu_lower
gamma_adj = adjust_shadowprice(gamma, 'Link', n)
revenue = (gamma_adj * lt.p0).sum()
assert close(investment, revenue)


# Storage Units

s = n.storage_units
st = n.storage_units_t


# This must be equal
assert close(s.capital_cost, (st.mu_upper_p_dispatch + st.mu_upper_p_store +
                              st.mu_upper_state_of_charge * s.max_hours).sum()
                             + s.mu_upper_p_nom + s.mu_lower_p_nom)

store, bus = s.index[0], s.iloc[0].bus
assert close((s.marginal_cost + st.mu_lower_p_dispatch + st.mu_upper_p_dispatch +
              st.mu_state_of_charge / s.efficiency_dispatch)[store],
             n.buses_t.marginal_price[bus])

investment = s.eval('p_nom_opt * capital_cost')

gamma = (st.mu_lower_p_dispatch + st.mu_upper_p_dispatch +
         st.mu_state_of_charge / s.efficiency_dispatch)
gamma_adj = adjust_shadowprice(gamma, 'StorageUnit', n)
revenue = (gamma_adj * st.p_dispatch - gamma * st.p_store).sum()
assert close(investment.where(st.p_dispatch.sum()>0, 0), revenue)




