# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 10:36:21 2019

@author: KatharinaG
"""
import argparse
import glob
import os
import xarray as xr

import sys
sys.path.append('../')

from merra2download import download_month

from dask.diagnostics import ProgressBar
ProgressBar().register()

from paths_usa import mer_path

parser = argparse.ArgumentParser(description='Insert date')
parser.add_argument('-date')
args = parser.parse_args()
ym = args.date

region = 'USA'

user = 'RE_EXTREME'
password = 'Re_extreme666!'
var = ['DISPH', 'U10M', 'U50M', 'V10M', 'V50M']

lat1 = 12
lat2 = 68
lon1 = -173
lon2 = -64

opath = mer_path

ofile = opath + '/merra2_wind_' + region + '_' + ym + '.nc'
if ofile in glob.glob(opath + '/*'):
    print(ofile, ' already there')
    exit()

if not os.path.exists(opath):
    print('Directory does not exist')
    exit()
elif not os.path.exists(opath + '/temp'):
    os.mkdir(opath + '/temp')

download_month(ym, lon1, lat1, lon2, lat2, var, user, password, opath + '/temp')

files = glob.glob(opath + '/temp/MERRA2_???.tavg1_2d_slv_Nx.' + ym + '??.nc4.nc')

d = xr.open_mfdataset(files)

print('merging month ' + ym + ' ...')
d.to_netcdf(ofile)

for file in files:
    os.remove(file)
