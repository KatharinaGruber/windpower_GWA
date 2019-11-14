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

print('correlating year ',year)

    
# calculate correlations for one year
if(results_path + '/cor_mer_'+year+'.nc' not in glob.glob(results_path + '/cor_mer_*.nc')):
    d_m = xr.open_dataset(results_path + '/oneyear/windpower_stat_MERRA2_'+year+'.nc')
    cor_m = np.corrcoef(d_m.wp,rowvar=False)
    xr.DataArray(cor_m,name='c').to_netcdf(results_path + '/cor_mer_'+year+'.nc')
    del(cor_m)

if(results_path + '/cor_era_'+year+'.nc' not in glob.glob(results_path + '/cor_era_*.nc')):
    d_e = xr.open_dataset(results_path + '/oneyear/windpower_stat_ERA5_'+year+'.nc')
    cor_e = np.corrcoef(d_e.wp,rowvar=False)
    xr.DataArray(cor_e,name='c').to_netcdf(results_path + '/cor_era_'+year+'.nc')
    del(cor_e)

if(results_path + '/cor_mer_gwa_'+year+'.nc' not in glob.glob(results_path + '/cor_mer_gwa_*.nc')):
    d_mg = xr.open_mfdataset(results_path + '/oneyear/windpower_*_MERRA2_GWA_'+year+'.nc')
    cor_mg = np.corrcoef(d_mg.wp,rowvar=False)
    xr.DataArray(cor_mg,name='c').to_netcdf(results_path + '/cor_mer_gwa_'+year+'.nc')
    del(cor_mg)

if(results_path + '/cor_era_gwa_'+year+'.nc' not in glob.glob(results_path + '/cor_era_gwa_*.nc')):
    d_eg = xr.open_mfdataset(results_path + '/oneyear/windpower_*_ERA5_GWA_'+year+'.nc')
    cor_eg = np.corrcoef(d_eg.wp,rowvar=False)
    xr.DataArray(cor_eg,name='c').to_netcdf(results_path + '/cor_era_gwa_'+year+'.nc')
    del(cor_eg)
