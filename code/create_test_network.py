#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 13:46:12 2020

@author: fabian
"""

import pypsa
import shutil

n = pypsa.Network(snakemake.input.network)

n.set_snapshots(n.snapshots[:50])

n.export_to_netcdf(snakemake.output.network)


shutil.copy(snakemake.input.regions, snakemake.output.regions)