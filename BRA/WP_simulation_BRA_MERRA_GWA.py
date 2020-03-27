import argparse
import datetime
import glob
import math
import numpy as np
import os
import pandas as pd
import rasterio
import statsmodels.api as sm
import time
import xarray as xr

import sys
sys.path.append('../')

from functools import reduce
from scipy.interpolate import interp1d
from utils import windpower_simulation_merra2

from dask.diagnostics import ProgressBar
ProgressBar().register()

from paths_bra import *

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

# Simulate wind power with MERRA-2
wind = xr.open_mfdataset(mer_path + "/eff_ws/merra2_wind_BRA_*.nc", chunks = {'time': 100})
alpha = xr.open_mfdataset(mer_path + "/eff_ws/merra2_alpha_BRA_*.nc", chunks = {'time': 100})

# with GWA
outfile = results_path + '/windpower_??_MERRA2_GWA.nc'
turbine_data = pd.read_csv(bra_path + '/turbine_data_mer_gwa' + GWA + '.csv', parse_dates=['commissioning'], index_col = 0)

if results_path + '/windpower_' + state + '_MERRA2_GWA.nc' not in glob.glob(outfile):
    print('calculating MERRA2 ' + state + ' GWA')
    if GWA == "3":
        GWA = xr.open_rasterio(bra_path+'/GWA/GWA3_BRA50m.tif')
    else:
        GWA = xr.open_dataarray(bra_path+'/GWA/GWA2_BRA50m.nc')
    ind = turbine_data.state == state
    wps = windpower_simulation_merra2(wind.wh50,
                                      alpha.alpha,
                                      turbine_data.height.values[ind],
                                      turbine_data.capacity.values[ind],
                                      turbine_data.sp.values[ind],
                                      turbine_data.lon.values[ind],
                                      turbine_data.lat.values[ind],
                                      turbine_data.commissioning.values[ind],
                                      GWA)
    # save as netcdf
    wps.drop(['x','y']).to_dataset(name='wp').to_netcdf(results_path+"/windpower_"+state+"_MERRA2_GWA.nc")
