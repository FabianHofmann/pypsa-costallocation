#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 11:06:12 2020

@author: fabian
"""

import pypsa
from helpers import adjust_shadowprice
import numpy as np

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('check_shadowprices', nname='de50bf')

if not 'test' in snakemake.input[0]:
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
    subsidy = g.p_nom_opt * g.mu_lower_p_nom
    scarcity = g.p_nom_max.replace(np.inf, 0) * g.mu_upper_p_nom
    assert close(investment, revenue + subsidy + scarcity)



    # Lines/Links

    l = n.lines
    lt = n.lines_t

    investment = l.eval('s_nom_opt * capital_cost')
    revenue = ((lt.mu_upper + lt.mu_lower) * lt.p0).sum()
    subsidy = l.s_nom_opt * l.mu_lower_s_nom
    scarcity = l.s_nom_opt * l.mu_upper_s_nom
    assert close(investment, revenue + subsidy + scarcity)


    l = n.links
    lt = n.links_t

    investment = l.eval('p_nom_opt * capital_cost')
    revenue = ((lt.mu_upper + lt.mu_lower) * lt.p0).sum()
    subsidy = l.p_nom_opt * l.mu_lower_p_nom
    scarcity = l.p_nom_opt * l.mu_upper_p_nom
    assert close(investment, revenue + subsidy + scarcity)



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


    gamma = (st.mu_lower_p_dispatch + st.mu_upper_p_dispatch +
             st.mu_state_of_charge / s.efficiency_dispatch)
    investment = s.eval('p_nom_opt * capital_cost')
    # gamma + marginal_cost is the lmp
    revenue = ((gamma) * st.p_dispatch - (gamma + s.marginal_cost) * st.p_store).sum()
    subsidy = s.p_nom_opt * s.mu_lower_p_nom
    scarcity = s.p_nom_opt * s.mu_upper_p_nom
    assert close(investment, revenue + subsidy + scarcity)





