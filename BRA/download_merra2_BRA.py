# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 10:36:21 2019

@author: KatharinaG
"""
import argparse
import glob
import numpy as np
import os
import xarray as xr

import sys
sys.path.append('../')

from merra2download import download_month

from dask.diagnostics import ProgressBar
ProgressBar().register()
from multiprocessing import Pool

from paths_bra import mer_path
from merra_cred import user, password


yms = np.char.replace(np.arange(np.datetime64('1987-01'),
                                np.datetime64('2020-01')).astype('str'),'-','')
region = 'BRA'

var = ['DISPH', 'U10M', 'U50M', 'V10M', 'V50M']

lat1 = -36
lat2 = 5.5
lon1 = -74.1
lon2 = -33

opath = mer_path



def dl_month_usa(ym):
    ofile = opath + '/merra2_wind_' + region + '_' + ym + '.nc'
    if ofile in glob.glob(opath + '/*'):
        print(ofile, ' already there')
        return()

    if not os.path.exists(opath):
        print('Directory does not exist')
        return()
    elif not os.path.exists(opath + '/temp'):
        os.mkdir(opath + '/temp')
    download_month(ym, lon1, lat1, lon2, lat2, var, user, password, opath + '/temp')

    files = glob.glob(opath + '/temp/MERRA2_???.tavg1_2d_slv_Nx.' + ym + '??.nc4.nc')

    d = xr.open_mfdataset(files)

    print('merging month ' + ym + ' ...')
    d.to_netcdf(ofile)

    for file in files:
        os.remove(file)
        
if __name__ == '__main__':
    pool = Pool()
    pool.map(dl_month_usa,yms)
