#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:24:39 2020

@author: fabian
"""

import pypsa
import pandas as pd

def n1_t1_g2_wo(d = 70, marginal_cost=[5, 10], p_nom=40):
    """
    $G_1 < d < G_1 + G_2$

    Note, for such a system the total system cost $TC$ are less then the nodal
    payments as soon as one generator is at its limit:

    $\lambda = o_s + \bar{\mu_s} \;\;\;\;\; \forall s$

    $d \, \lambda = d \, (o_s + \bar{\mu_s} ) \ge \sum_s g_s \, o_s = TC$
    """
    n = pypsa.Network()
    n.add('Bus', 'Bus1')
    n.add('Load', 0, bus='Bus1', p_set=d)
    n.madd('Generator', ['Gen0', 'Gen1'], bus='Bus1',
           marginal_cost=marginal_cost, p_nom=p_nom)
    return n



def n1_t1_g2_w(d = 70, marginal_cost=[5, 10], capital_cost=[100, 30]):
    """
    Here the total cost are payed by the consumers

    $TC = \lambda d$
    """
    n = pypsa.Network()
    n.add('Bus', 'Bus1')
    n.add('Load', 0, bus='Bus1', p_set=d)
    n.madd('Generator', ['Gen0', 'Gen1'], bus='Bus1', p_nom_extendable=True,
           marginal_cost=marginal_cost, capital_cost=capital_cost)
    return n


def n1_t9_g2_w(d = 70, marginal_cost=[2, 4.], capital_cost=[50, 30]):
    """Initialize network with base and one peak demand at one bus."""
    load = pd.DataFrame({0: [50, 49, 48, 47, 46, 90, 46.5, 47.5, 48.5]})

    n = pypsa.Network()
    n.set_snapshots(list(range(len(load))))
    n.add('Bus', 'Bus1')
    n.madd('Load', [0], bus='Bus1', p_set=load)
    n.madd('Generator', ['Gen0', 'Gen1'], bus='Bus1', p_nom_extendable=True,
           marginal_cost=marginal_cost, capital_cost=capital_cost)
    return n


def n2_t9_g2_w(d = 70, marginal_cost=[1, 2.4], capital_cost=[10, 4]):
    """Initialize network with base and one peak demand at one bus."""
    load = pd.DataFrame({0: [10, 10, 10, 11]})

    n = pypsa.Network()
    n.set_snapshots(list(range(len(load))))
    n.add('Bus', 'Bus1')
    n.add('Bus', 'Bus2')
    n.madd('Load', [0], bus='Bus1', p_set=load)
    n.madd('Generator', ['Base', 'Peak'], bus=['Bus1', 'Bus2'],
           p_nom_extendable=True,
           marginal_cost=marginal_cost, capital_cost=capital_cost)
    n.madd('Line', ['1'], s_nom_extendable=True, x=0.01, bus0='Bus1',
           bus1='Bus2', capital_cost=1)
    return n


def n2_t1_g2_w(n0='1', n1='2', d0=60, d1=90, o0=50, o1=200, c0=500, c1=500,
               c_l=100):
    n = pypsa.Network()

    n.madd('Bus', [n0, n1], x=[0, 1], y=[0, 0.4])
    n.madd('Load', [0, 1], bus=[n0, n1], p_set=[d0, d1])
    n.madd('Generator', [0, 1], bus=[n0, n1], p_nom_extendable=True,
           marginal_cost=[o0, o1], capital_cost=[c0, c1], p_nom_max=100)
    n.madd('Line', ['1'], s_nom_extendable=True, x=0.01, bus0=[n0], bus1=[n1],
           capital_cost=c_l)
    return n


def n_ac_dc():
    import netallocation as ntl

    n = ntl.test.get_network_ac_dc()
    n.generators['p_nom'] = 0
    n.generators['p_nom_min'] = 0
    n.links['p_nom'] = 0
    n.lines['s_nom'] = 0
    n.set_snapshots(n.snapshots[[0]])
    return n