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

from paths_bra import *


# Simulate wind power with ERA5
wind = xr.open_mfdataset(era_path + "/eff_ws/era5_wind_BRA_*.nc", chunks = {'time': 100})
alpha = xr.open_mfdataset(era_path + "/eff_ws/era5_alpha_BRA_*.nc", chunks = {'time': 100})

# without GWA
outfile = results_path + '/windpower_stat_ERA5.nc'
turbine_data = pd.read_csv(bra_path + '/turbine_data_era.csv', parse_dates=['commissioning'], index_col = 0)

if outfile not in glob.glob(results_path+'/*'):
	wps = windpower_simulation_era5(wind.wh100,
                                    alpha.alpha,
                                    turbine_data.height.values,
                                    turbine_data.capacity.values,
                                    turbine_data.sp.values,
                                    turbine_data.lon.values,
                                    turbine_data.lat.values,
                                    turbine_data.commissioning.values)
	# save as netcdf
	wps.to_dataset(name='wp').to_netcdf(outfile)