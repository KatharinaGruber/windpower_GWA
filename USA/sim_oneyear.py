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
from utils import windpower_simulation_era5_large
from utils import windpower_simulation_merra2_large


from dask.diagnostics import ProgressBar
ProgressBar().register()

from paths_usa import *

parser = argparse.ArgumentParser()
parser.add_argument('-year',help='Insert year')
args = parser.parse_args()
year = args.year

print('calculating year ',year)

# prepare turbine data
if len(glob.glob(usa_path+'/turbine_data_*.csv'))!=4:
    exec(open('prepare_USA_turbines.py').read())

# prepare reanalysis data: calcualte effective wind speeds and alpha
if len(glob.glob(mer_path + '/eff_ws/*')) != 4:
    exec(open('prepare_USA_MERRA2.py').read())
if len(glob.glob(era_path + '/eff_ws/*')) != 18:
    exec(open('prepare_USA_ERA5.py').read())


# Get data
wind_merra = xr.open_mfdataset(mer_path + "/eff_ws/merra2_wind_USA_*.nc", chunks = {'time': 38}).sel(time=np.arange(year+'-01-01',str(int(year)+1)+'-01-01',dtype='datetime64[h]'))
alpha_merra = xr.open_mfdataset(mer_path + "/eff_ws/merra2_alpha_USA_*.nc", chunks = {'time': 38}).sel(time=np.arange(year+'-01-01',str(int(year)+1)+'-01-01',dtype='datetime64[h]'))
wind_era = xr.open_mfdataset(era_path + "/eff_ws/era5_wind_USA_*.nc", chunks = {'time': 38}).sel(time=np.arange(year+'-01-01',str(int(year)+1)+'-01-01',dtype='datetime64[h]'))
alpha_era = xr.open_mfdataset(era_path + "/eff_ws/era5_alpha_USA_*.nc", chunks = {'time': 38}).sel(time=np.arange(year+'-01-01',str(int(year)+1)+'-01-01',dtype='datetime64[h]'))


turbine_data_mer = pd.read_csv(usa_path + '/turbine_data_mer.csv', parse_dates=['commissioning'])
turbine_data_era = pd.read_csv(usa_path + '/turbine_data_era.csv', parse_dates=['commissioning'])
turbine_data_mer_gwa = pd.read_csv(usa_path + '/turbine_data_mer_gwa.csv', parse_dates=['commissioning'])
turbine_data_era_gwa = pd.read_csv(usa_path + '/turbine_data_era_gwa.csv', parse_dates=['commissioning'])

# without GWA
outfile = results_path + '/oneyear/windpower_stat_MERRA2_'+year+'.nc'
if outfile not in glob.glob(results_path+'/oneyear/*'):
    wps = windpower_simulation_merra2(wind_merra.wh50,
                                      alpha_merra.alpha,
                                      turbine_data_mer.height.values,
                                      turbine_data_mer.capacity.values,
                                      turbine_data_mer.sp.values,
                                      turbine_data_mer.lon.values,
                                      turbine_data_mer.lat.values,
                                      pd.to_datetime([year+'-01-01']*len(turbine_data_mer.commissioning)))
    # save as netcdf
    wps.to_dataset(name='wp').to_netcdf(outfile)
    del(wps)
    
outfile = results_path + '/oneyear/windpower_stat_ERA5_'+year+'.nc'
if outfile not in glob.glob(results_path+'/oneyear/*'):
    wps = windpower_simulation_era5(wind_era.wh100,
                                    alpha_era.alpha,
                                    turbine_data_era.height.values,
                                    turbine_data_era.capacity.values,
                                    turbine_data_era.sp.values,
                                    turbine_data_era.lon.values,
                                    turbine_data_era.lat.values,
                                    pd.to_datetime([year+'-01-01']*len(turbine_data_era.commissioning)))
    # save as netcdf
    wps.to_dataset(name='wp').to_netcdf(outfile)
    del(wps)



os.chdir(script_path)
os.system('bash run_USA_MERRA_GWA_simulation_oneyear.sh 8 '+year)
os.system('bash run_USA_ERA_GWA_simulation_oneyear.sh 8 '+year)
