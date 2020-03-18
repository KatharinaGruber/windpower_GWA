# more or less the same simulation, but split up in to chunks that fit into memory
# for large states (CA, IA, KS, OK, TX)
# chunks of size 2000 (2000 locations in one part calculated and create temporary file)
location_chunk = 2000

import argparse
import datetime
import glob
import math
import numpy as np
import os
import rasterio
import statsmodels.api as sm
import pandas as pd
import statsmodels.api as sm
import time
import xarray as xr

import sys
sys.path.append('../')

from functools import reduce
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from utils import windpower_simulation_era5_large

from dask.diagnostics import ProgressBar
ProgressBar().register()

from paths_usa import *

# get state and GWA version
parser = argparse.ArgumentParser(description='Insert state and optionally GWA')
parser.add_argument('-state')
parser.add_argument('-GWA')
args = parser.parse_args()
state = args.state
if(args.GWA == None):
    GWA = 3
else:
    GWA = args.GWA

if GWA == "2":
    results_path = results_path + '/results_GWA2'


outfile = results_path + '/windpower_??_ERA5_GWA.nc'
if results_path + '/windpower_' + state + '_ERA5_GWA.nc' not in glob.glob(outfile):

    wind = xr.open_mfdataset(era_path + "/eff_ws/era5_wind_USA_*.nc", chunks = {'time': 38})
    alpha = xr.open_mfdataset(era_path + "/eff_ws/era5_alpha_USA_*.nc", chunks = {'time': 38})
    # with GWA
    turbine_data_era_gwa = pd.read_csv(usa_path + '/turbine_data_era_gwa.csv', parse_dates=['commissioning'])
    if GWA == 3:
        if state == 'PR':
            GWA = xr.open_rasterio(usa_path+'/GWA/GWA3_PR100m.tif')
        else:
            GWA = xr.open_rasterio(usa_path+'/GWA/GWA3_USA100m.tif')
    else:
        if state == 'AK':
            GWA = xr.open_rasterio(usa_path+'/GWA/GWA_AK100m.tif')
        elif state == 'HI':
            GWA = xr.open_rasterio(usa_path+'/GWA/GWA_HI100m.tif')
        elif state == 'PR':
            GWA = xr.open_rasterio(usa_path+'/GWA/GWA_PR100m.tif')
        else:
            GWA = xr.open_rasterio(usa_path+'/GWA/GWA_USA100m.tif')
    ind = turbine_data_era_gwa.state == state

    print('calculating ERA5 ' + state + ' GWA')

    t1 = time.time()
    # number of locations in state
    dat_len = sum(ind)
    numit = round(dat_len/location_chunk+0.5) # number of necessary iterations
    i1 = 0
    i2 = i1 + location_chunk
    for it in range(numit):
        outfile_temp = results_path + "/wp_"+state+"_ERA5_GWA_temp" + str(it+1) +".nc"
        if i2 > dat_len:
            i2 = dat_len
        if outfile_temp not in glob.glob(results_path + "/wp_"+state+"_ERA5_GWA_temp*.nc"):
            wps = windpower_simulation_era5_large(wind.wh100,
                                                  alpha.alpha,
                                                  turbine_data_era_gwa.height[ind].values[i1:i2],
                                                  turbine_data_era_gwa.capacity[ind].values[i1:i2],
                                                  turbine_data_era_gwa.sp[ind].values[i1:i2],
                                                  turbine_data_era_gwa.lon[ind].values[i1:i2],
                                                  turbine_data_era_gwa.lat[ind].values[i1:i2],
                                                  pd.to_datetime(turbine_data_era_gwa.commissioning[ind].values[i1:i2]).year.values,
                                                  GWA)
            # adapt numbers of locations in dataset
            wps = wps.assign_coords(location = np.arange(i1,i2))
            # save temporary file
            print('saving to '+results_path + "/wp_"+state+"_ERA5_GWA_temp" + str(it+1) +".nc")
            
            wps.to_dataset(name='wp').to_netcdf(results_path + "/wp_"+state+"_ERA5_GWA_temp" + str(it+1) +".nc")
            
            
            print('saved to '+results_path + "/wp_"+state+"_ERA5_GWA_temp" + str(it+1) +".nc")
            
            wps.close()
            del(wps)
        i1 = i2
        i2 = i2 + location_chunk
        print(round(i1/dat_len,3)*100,'% done in ',state)
        
    # merge  and delete temporary files
    wps = xr.open_mfdataset(results_path + "/wp_"+state+"_ERA5_GWA_temp*.nc", chunks = {'time': 100})
    
    print('saving to'+results_path + "/windpower_"+state+"_ERA5_GWA.nc")
    
    wps.drop(['x','y']).to_netcdf(results_path + "/windpower_"+state+"_ERA5_GWA.nc")
    
    print('saved to '+results_path + "/windpower_"+state+"_ERA5_GWA.nc")
    
    t2 = time.time()
    
    # remove temporary files
    for file in glob.glob(results_path + "/wp_"+state+"_ERA5_GWA_temp*.nc"):
        os.remove(file)


    print(t2-t1)



