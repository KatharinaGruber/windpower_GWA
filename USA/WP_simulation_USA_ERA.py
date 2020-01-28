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
from utils import windpower_simulation_era5_large


from dask.diagnostics import ProgressBar
ProgressBar().register()

from paths_usa import *


# prepare turbine data
if len(glob.glob(usa_path+'/turbine_data_*.csv'))!=4:
	exec(open('prepare_USA_turbines.py').read())
	
# prepare reanalysis data: calcualte effective wind speeds and alpha
if len(glob.glob(era_path + '/eff_ws/*')) != 18:
	exec(open('prepare_USA_ERA5.py').read())


# Simulate wind power with ERA5
wind = xr.open_mfdataset(era_path + "/eff_ws/era5_wind_USA_*.nc", chunks = {'time': 38})
alpha = xr.open_mfdataset(era_path + "/eff_ws/era5_alpha_USA_*.nc", chunks = {'time': 38})

# without GWA
outfile = results_path + '/windpower_stat_ERA5.nc'
turbine_data_era = pd.read_csv(usa_path + '/turbine_data_era.csv', parse_dates=['commissioning'])
if outfile not in glob.glob(results_path+'/*'):
	wps = windpower_simulation_era5_large(wind.wh100,
                                          alpha.alpha,
                                          turbine_data_era.height.values,
                                          turbine_data_era.capacity.values,
                                          turbine_data_era.sp.values,
                                          turbine_data_era.lon.values,
                                          turbine_data_era.lat.values,
                                          pd.to_datetime(turbine_data_era.commissioning.values).year.values)
	# save as netcdf
	wps.to_dataset(name='wp').to_netcdf(outfile)