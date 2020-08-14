#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 11:06:12 2020

@author: fabian
"""

import pypsa


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from _helpers import mock_snakemake
        snakemake = mock_snakemake('check_shadowprices', nname='de10gf')


n = pypsa.Network(snakemake.input[0])

tol = (n.objective + n.objective_constant)*1e-8


close = lambda df1, df2: ((df1 - df2).abs() < tol ).all()

# Generators

g = n.generators
gt = n.generators_t

assert close(g.capital_cost, (gt.mu_upper.mul(gt.p_max_pu, fill_value=1)).sum()
                              + g.mu_upper_p_nom + g.mu_lower_p_nom)

investment = g.eval('p_nom_opt * capital_cost')
revenue = (gt.mu_upper * gt.p).sum()

correction_factor_upper = g.capital_cost / (g.capital_cost - g.mu_upper_p_nom)
correction_term_lower = g.mu_lower_p_nom * g.p_nom_opt / gt.p.sum()

revenue_corrected = ((gt.mu_upper * correction_factor_upper
                      + correction_term_lower) * gt.p).sum()

assert close(investment.where(gt.p.sum()>0, 0), revenue_corrected)


# Lines/Links




# Storage Units

s = n.storage_units
st = n.storage_units_t


# This must be equal
assert close(s.capital_cost, (st.mu_upper_p_dispatch + st.mu_upper_p_store +
                              st.mu_upper_state_of_charge * s.max_hours).sum())

store, bus = s.index[0], s.iloc[0].bus
assert close((s.marginal_cost + st.mu_lower_p_dispatch + st.mu_upper_p_dispatch + st.mu_state_of_charge / s.efficiency_dispatch)[store], n.buses_t.marginal_price[bus])


