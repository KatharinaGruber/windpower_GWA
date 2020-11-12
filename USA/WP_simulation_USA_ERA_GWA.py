import argparse
import datetime
import glob
import math
import numpy as np
import os
import pandas as pd
import rasterio
import seaborn as sns
import statsmodels.api as sm
import time
import xarray as xr

import sys
sys.path.append('../')

from functools import reduce
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

from utils import power_curve
from utils import windpower_simulation_era5
from utils import windpower_simulation_merra2


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
    GWA = "3"
else:
    GWA = args.GWA

if GWA == "2":
    results_path = results_path + '/results_GWA2'
    if not os.path.exists(results_path):
        os.mkdir(results_path)
    startGWA = '1987'
    endGWA = '2016'
else:
    startGWA = '2008'
    endGWA = '2017'
# define start date for simulation
startyearmonth = '2000-12'

# Simulate wind power with ERA5
wind = xr.open_mfdataset(era_path + "/eff_ws/era5_wind_USA_*.nc", chunks = {'time': 38})
alpha = xr.open_mfdataset(era_path + "/eff_ws/era5_alpha_USA_*.nc", chunks = {'time': 38})

# with GWA
outfile = results_path + '/windpower_??_ERA5_GWA.nc'
turbine_data_era_gwa = pd.read_csv(usa_path + '/turbine_data_era_gwa' + GWA + '.csv', parse_dates=['commissioning'])
if results_path + '/windpower_' + state + '_ERA5_GWA.nc' not in glob.glob(outfile):
    print('calculating ERA5 ' + state + ' GWA')
    if GWA == "3":
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
    wps = windpower_simulation_era5(wind.wh100,
                                    alpha.alpha,
                                    turbine_data_era_gwa.height[ind].values,
                                    turbine_data_era_gwa.capacity[ind].values,
                                    turbine_data_era_gwa.sp[ind].values,
                                    turbine_data_era_gwa.lon[ind].values,
                                    turbine_data_era_gwa.lat[ind].values,
                                    pd.to_datetime(turbine_data_era_gwa.commissioning[ind].values).year.values,
                                    startyearmonth,
                                    GWA,
                                    startGWA,
                                    endGWA)
    # save as netcdf
    wps.drop(['x','y']).to_dataset(name='wp').to_netcdf(results_path+"/windpower_"+state+"_ERA5_GWA.nc")
