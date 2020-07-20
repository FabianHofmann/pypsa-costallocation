#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 13:46:12 2020

@author: fabian
"""

import pypsa

n = pypsa.Network(snakemake.input[0])

n.set_snapshots(n.snapshots[:10])

n.export_to_netcdf(snakemake.output[0])