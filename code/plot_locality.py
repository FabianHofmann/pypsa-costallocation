#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 17:47:01 2020

@author: fabian
"""
import netallocation as ntl
import pypsa
import xarray as xr
import geopandas as gpd
import cartopy.crs as ccrs
from config import to_explanation
import matplotlib.pyplot as plt
import os
from helpers import combine_oneports, scfmt, load

if 'snakemake' not in globals():
    from _helpers import mock_snakemake
    snakemake = mock_snakemake('plot_sparcity_maps', nname='test-de10bf',
                               method='ptpf', power='net')


n = pypsa.Network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.regions).set_index('name')
cost = (combine_oneports(xr.open_dataset(snakemake.input.costs).sum('snapshot'))
        .set_index({'branch': ['component', 'branch_i']})
        .transpose('sink', 'source', 'source_carrier', 'branch'))
