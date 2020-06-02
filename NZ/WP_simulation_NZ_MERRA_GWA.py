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

from paths_nz import *

# get GWA version
parser = argparse.ArgumentParser(description='Insert optionally GWA')
parser.add_argument('-GWA')
args = parser.parse_args()
if(args.GWA == None):
    GWA = "3"
else:
    GWA = args.GWA

if GWA == "2":
    results_path = results_path + '/results_GWA2'
    if not os.path.exists(results_path):
        os.mkdir(results_path)

# Simulate wind power with MERRA-2
wind = xr.open_mfdataset(mer_path + "/eff_ws/merra2_wind_NZ_*.nc", chunks = {'time': 46}).sel(time=slice('1997','2020'))
alpha = xr.open_mfdataset(mer_path + "/eff_ws/merra2_alpha_NZ_*.nc", chunks = {'time': 46}).sel(time=slice('1997','2020'))

# load windpark data
windparks = pd.read_csv(nz_path + "/windparks_NZ.csv", delimiter=';', parse_dates=['commissioning'])
# calculate specific power of turbines (in W)
windparks['sp'] = windparks.turb_cap*10**6/(windparks.d_rotor**2*np.pi/4)

# with GWA
outfile = results_path + '/windpower_??_MERRA2_GWA.nc'

if results_path + '/windpower_NZ_MERRA2_GWA.nc' not in glob.glob(outfile):
    print('calculating MERRA2 NZ GWA')
    if GWA == "3":
        GWA = xr.open_rasterio(nz_path+'/GWA/GWA3_NZ50m.tif')
    else:
        GWA = xr.open_dataarray(nz_path+'/GWA/GWA2_NZ50m.nc')
    wps = windpower_simulation_merra2(wind.wh50,
                                      alpha.alpha,
                                      windparks.Height.values,
                                      windparks.Capacity.values,
                                      windparks.sp.values,
                                      windparks.Longitude.values,
                                      windparks.Latitude.values,
                                      windparks.commissioning.values,
                                      GWA)

    # save as netcdf
    wps.drop(['x','y']).to_dataset(name='wp').to_netcdf(results_path+"/windpower_NZ_MERRA2_GWA.nc")
