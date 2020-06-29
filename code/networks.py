#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:24:39 2020

@author: fabian
"""

import pypsa
import pandas as pd

def n1_t1_g2_wo(d = 70, marginal_cost=[5, 10], p_nom=40):
    n = pypsa.Network()
    n.add('Bus', 'Bus1')
    n.add('Load', 0, bus='Bus1', p_set=d)
    n.madd('Generator', ['Gen0', 'Gen1'], bus='Bus1',
           marginal_cost=marginal_cost, p_nom=p_nom)
    return n



def n1_t1_g2_w(d = 70, marginal_cost=[5, 10], capital_cost=[100, 30]):
    n = pypsa.Network()
    n.add('Bus', 'Bus1')
    n.add('Load', 0, bus='Bus1', p_set=d)
    n.madd('Generator', ['Gen0', 'Gen1'], bus='Bus1', p_nom_extendable=True,
           marginal_cost=marginal_cost, capital_cost=capital_cost)
    return n


def n1_t9_g2_w(d = 70, marginal_cost=[5, 10], capital_cost=[100, 30]):
    """Initialize network with base and one peak demand at one bus."""
    load = pd.DataFrame([[50, 50, 50, 50, 50, 90, 50, 50, 50]], index=[0])

    n = pypsa.Network()
    n.set_snapshots(range(len(load)))
    n.add('Bus', 'Bus1')
    n.madd('Load', [0], bus='Bus1', p_set=load)
    n.madd('Generator', ['Gen0', 'Gen1'], bus='Bus1', p_nom_extendable=True,
           marginal_cost=marginal_cost, capital_cost=capital_cost)
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
    n.links['p_nom'] = 0
    n.lines['s_nom'] = 0
    n.set_snapshots(n.snapshots[[0]])
    return n