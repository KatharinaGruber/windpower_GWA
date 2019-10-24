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

parser = argparse.ArgumentParser(description='Insert state')
parser.add_argument('-state')
args = parser.parse_args()
state = args.state


# prepare turbine data
if len(glob.glob(usa_path+'/turbine_data_*.csv'))!=4:
	exec(open('prepare_USA_turbines.py').read())
	
# prepare reanalysis data: calcualte effective wind speeds and alpha
if len(glob.glob(era_path + '/eff_ws/*')) != 18:
	exec(open('prepare_USA_ERA5.py').read())


# Simulate wind power with ERA5
wind = xr.open_mfdataset(era_path + "/eff_ws/era5_wind_USA_*.nc", chunks = {'time': 38})
alpha = xr.open_mfdataset(era_path + "/eff_ws/era5_alpha_USA_*.nc")

# with GWA
outfile = results_path + '/windpower_??_ERA5_GWA.nc'
turbine_data_era_gwa = pd.read_csv(usa_path + '/turbine_data_era_gwa.csv', parse_dates=['commissioning'])
if results_path + '/windpower_' + state + '_ERA5_GWA.nc' not in glob.glob(outfile):
	print('calculating ERA5 ' + state + ' GWA')
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
									turbine_data_era_gwa.commissioning[ind].values,
									GWA)
	# save as netcdf
	wps.to_dataset(name='wp').to_netcdf(results_path+"/windpower_"+state+"_ERA5_GWA.nc")
