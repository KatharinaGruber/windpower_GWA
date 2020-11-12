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
from utils import windpower_simulation_era5

from dask.diagnostics import ProgressBar
ProgressBar().register()

from paths_zaf import *

if not os.path.isdir(results_path):
    os.mkdir(results_path)

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
    startGWA = '1987'
    endGWA = '2016'
else:
    startGWA = '2008'
    endGWA = '2017'
# define start date for simulation
startyear = '2013'

# Simulate wind power with ERA5
wind = xr.open_mfdataset(era_path + "/eff_ws/era5_wind_ZAF_*.nc", chunks = {'time': 46})
alpha = xr.open_mfdataset(era_path + "/eff_ws/era5_alpha_ZAF_*.nc", chunks = {'time': 46})

# load windpark data
windparks = pd.read_csv(zaf_path + "/windparks_ZAF.csv", parse_dates=['commissioning'])
# calculate specific power of turbines (in W)
windparks['sp'] = windparks.Turb_cap*10**6/(windparks.Diameter**2*np.pi/4)

# with GWA
outfile = results_path + '/windpower_ZAF_ERA5_GWA.nc'
if outfile not in glob.glob(results_path + '/*.nc'):
    print('calculating ERA5 ZAF GWA')
    if GWA == "3":
        GWA = xr.open_rasterio(zaf_path+'/GWA/GWA3_ZAF100m.tif')
    else:
        GWA = xr.open_dataarray(zaf_path+'/GWA/GWA2_ZAF100m.nc')
    wps = windpower_simulation_era5(wind.wh100,
                                    alpha.alpha,
                                    windparks.Height.values,
                                    windparks.Capacity.values,
                                    windparks.sp.values,
                                    windparks.Longitude.values,
                                    windparks.Latitude.values,
                                    windparks.commissioning.values,
                                    startyear,
                                    GWA,
                                    startGWA,
                                    endGWA)
    # save as netcdf
    wps.drop(['x','y']).to_dataset(name='wp').to_netcdf(results_path+"/windpower_ZAF_ERA5_GWA.nc")
