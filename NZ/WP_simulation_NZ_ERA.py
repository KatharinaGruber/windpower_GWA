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

from utils import power_curve
from utils import windpower_simulation_era5


from dask.diagnostics import ProgressBar
ProgressBar().register()

from paths_nz import *


# Simulate wind power with ERA5
wind = xr.open_mfdataset(era_path + "/eff_ws/era5_wind_NZ_*.nc", chunks = {'time': 100}).sel(time=slice('1997','2020'))
alpha = xr.open_mfdataset(era_path + "/eff_ws/era5_alpha_NZ_*.nc", chunks = {'time': 100}).sel(time=slice('1997','2020'))

# load windpark data
windparks = pd.read_csv(nz_path + "/windparks_NZ.csv", delimiter=';', parse_dates=['commissioning'])
# calculate specific power of turbines (in W)
windparks['sp'] = windparks.turb_cap*10**6/(windparks.d_rotor**2*np.pi/4)

# without GWA
outfile = results_path + '/windpower_NZ_ERA5.nc'

if outfile not in glob.glob(results_path+'/*'):
	wps = windpower_simulation_era5(wind.wh100,
                                    alpha.alpha,
                                    windparks.Height.values,
                                    windparks.Capacity.values,
                                    windparks.sp.values,
                                    windparks.Longitude.values,
                                    windparks.Latitude.values,
                                    windparks.commissioning.values)
	# save as netcdf
	wps.to_dataset(name='wp').to_netcdf(outfile)