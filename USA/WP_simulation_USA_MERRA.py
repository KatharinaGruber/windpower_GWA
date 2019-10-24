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

# prepare turbine data
if len(glob.glob(usa_path+'/turbine_data_*.csv'))!=4:
	exec(open('prepare_USA_turbines.py').read())
	
# prepare reanalysis data: calcualte effective wind speeds and alpha
if len(glob.glob(mer_path + '/eff_ws/*')) != 4:
	exec(open('prepare_USA_MERRA2.py').read())
if len(glob.glob(era_path + '/eff_ws/*')) != 18:
	exec(open('prepare_USA_ERA5.py').read())


# Simulate wind power with MERRA-2
wind = xr.open_mfdataset(mer_path + "/eff_ws/merra2_wind_USA_*.nc", chunks = {'time': 38})
alpha = xr.open_mfdataset(mer_path + "/eff_ws/merra2_alpha_USA_*.nc")

# without GWA
outfile = results_path + '/windpower_stat_MERRA2.nc'
turbine_data_mer = pd.read_csv(usa_path + '/turbine_data_mer.csv', parse_dates=['commissioning'])
if outfile not in glob.glob(results_path+'/*'):
	wps = windpower_simulation_merra2(wind.wh50,
									  alpha.alpha,
									  turbine_data_mer.height.values,
									  turbine_data_mer.capacity.values,
									  turbine_data_mer.sp.values,
									  turbine_data_mer.lon.values,
									  turbine_data_mer.lat.values,
									  turbine_data_mer.commissioning.values)
	# save as netcdf
	wps.to_dataset(name='wp').to_netcdf(outfile)