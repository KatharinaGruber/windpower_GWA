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


# Simulate wind power with MERRA-2
wind = xr.open_mfdataset(mer_path + "/eff_ws/merra2_wind_USA_*.nc", chunks = {'time': 38})
alpha = xr.open_mfdataset(mer_path + "/eff_ws/merra2_alpha_USA_*.nc", chunks = {'time': 38})

# with GWA
outfile = results_path + '/windpower_??_MERRA2_GWA.nc'
#states = np.unique(turbine_data_mer.state)
turbine_data_mer_gwa = pd.read_csv(usa_path + '/turbine_data_mer_gwa' + GWA + '.csv', parse_dates=['commissioning'])

if results_path + '/windpower_' + state + '_MERRA2_GWA.nc' not in glob.glob(outfile):
    print('calculating MERRA2 ' + state + ' GWA')
    if GWA == "3":
        if state == 'PR':
            GWA = xr.open_rasterio(usa_path+'/GWA/GWA3_PR50m.tif')
        else:
            GWA = xr.open_rasterio(usa_path+'/GWA/GWA3_USA50m.tif')
        ind = turbine_data_mer_gwa.state == state
    else:
        if state == 'AK':
            GWA = xr.open_rasterio(usa_path+'/GWA/GWA_AK50m.tif')
        elif state == 'HI':
            GWA = xr.open_rasterio(usa_path+'/GWA/GWA_HI50m.tif')
        elif state == 'PR':
            GWA = xr.open_rasterio(usa_path+'/GWA/GWA_PR50m.tif')
        else:
            GWA = xr.open_rasterio(usa_path+'/GWA/GWA_USA50m.tif')
    ind = turbine_data_mer_gwa.state == state
    wps = windpower_simulation_merra2(wind.wh50,
                                      alpha.alpha,
                                      turbine_data_mer_gwa.height[ind].values,
                                      turbine_data_mer_gwa.capacity[ind].values,
                                      turbine_data_mer_gwa.sp[ind].values,
                                      turbine_data_mer_gwa.lon[ind].values,
                                      turbine_data_mer_gwa.lat[ind].values,
                                      pd.to_datetime(turbine_data_mer_gwa.commissioning[ind].values).year.values,
                                      GWA)
    # save as netcdf
    wps.drop(['x','y']).to_dataset(name='wp').to_netcdf(results_path+"/windpower_"+state+"_MERRA2_GWA.nc")