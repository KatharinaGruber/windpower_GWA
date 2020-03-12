# Script for preparing MERRA-2 reanalysis data:
# 1. Calculate effective wind speeds from u and v wind speeds in two heights (10 and 50m)
# 2. Calculate alpha friction coefficient

import glob
import numpy as np
import xarray as xr
from paths_usa import mer_path

# first period: 2000-2009
if len(glob.glob(mer_path + "/eff_ws/merra2_*_USA_2000-2009.nc")) < 2:
    data = xr.open_mfdataset(mer_path + "/merra2_wind_USA_200*.nc", chunks = {'time': 38})
    wh10 = ((data.U10M**2+data.V10M**2)**0.5).compute()
    wh50 = ((data.U50M**2+data.V50M**2)**0.5).compute()

    eff_ws = xr.Dataset({'wh10': wh10,
                         'wh50': wh50})
                     
    eff_ws.to_netcdf(mer_path+"/eff_ws/merra2_wind_USA_2000-2009.nc")

    alpha = (xr.ufuncs.log(eff_ws.wh50/eff_ws.wh10)/np.log(50/(10+data.DISPH))).compute()

    xr.Dataset({'alpha': alpha}).to_netcdf(mer_path+"/eff_ws/merra2_alpha_USA_2000-2009.nc")

# second period: 2010-2019
if len(glob.glob(mer_path + "/eff_ws/merra2_*_USA_2010-2019.nc")) < 2:
    data = xr.open_mfdataset(mer_path + "/merra2_wind_USA_201*.nc", chunks = {'time': 38})
    wh10 = ((data.U10M**2+data.V10M**2)**0.5).compute()
    wh50 = ((data.U50M**2+data.V50M**2)**0.5).compute()

    eff_ws = xr.Dataset({'wh10': wh10,
                         'wh50': wh50})
                     
    eff_ws.to_netcdf(mer_path+"/eff_ws/merra2_wind_USA_2010-2019.nc")

    alpha = (xr.ufuncs.log(eff_ws.wh50/eff_ws.wh10)/np.log(50/(10+data.DISPH))).compute()

    xr.Dataset({'alpha': alpha}).to_netcdf(mer_path+"/eff_ws/merra2_alpha_USA_2010-2019.nc")