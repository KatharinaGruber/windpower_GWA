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
from utils import windpower_simulation_merra2


from dask.diagnostics import ProgressBar
ProgressBar().register()

from paths_zaf import *

if not os.path.isdir(results_path):
    os.mkdir(results_path)

# Simulate wind power with MERRA-2
wind = xr.open_mfdataset(mer_path + "/eff_ws/merra2_wind_ZAF_*.nc", chunks = {'time': 46})
alpha = xr.open_mfdataset(mer_path + "/eff_ws/merra2_alpha_ZAF_*.nc", chunks = {'time': 46})

# load windpark data
windparks = pd.read_csv(zaf_path + "/windparks_ZAF.csv", parse_dates=['commissioning'])
# calculate specific power of turbines (in W)
windparks['sp'] = windparks.Turb_cap*10**6/(windparks.Diameter**2*np.pi/4)

# without GWA
outfile = results_path + '/windpower_ZAF_MERRA2.nc'
if outfile not in glob.glob(results_path+'/*'):
	wps = windpower_simulation_merra2(wind.wh50,
                                      alpha.alpha,
                                      windparks.Height.values,
                                      windparks.Capacity.values,
                                      windparks.sp.values,
                                      windparks.Longitude.values,
                                      windparks.Latitude.values,
                                      windparks.commissioning.values)
	# save as netcdf
	wps.to_dataset(name='wp').to_netcdf(outfile)