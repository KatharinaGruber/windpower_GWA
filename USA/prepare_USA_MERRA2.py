# Script for preparing MERRA-2 reanalysis data:
# 1. Calculate effective wind speeds from u and v wind speeds in two heights (10 and 50m)
# 2. Calculate alpha friction coefficient

import glob
import numpy as np
import xarray as xr
from paths_usa import mer_path

for dec in ['198','199','200','201']:
    if dec=='198':
        year1 = '1987'
    else:
        year1 = dec + '0'
    year2 = dec + '9'
    if len(glob.glob(mer_path + '/eff_ws/merra2_*_USA_' + year1 + '-' + year2 + '.nc'))<2:
        print('calculating wind ' + year1 + '-' + year2)
        data = xr.open_mfdataset(mer_path + "/merra2_wind_USA_" + dec + "*.nc", chunks = {'time': 38})
        wh10 = ((data.U10M**2+data.V10M**2)**0.5).compute()
        wh50 = ((data.U50M**2+data.V50M**2)**0.5).compute()
        
        eff_ws = xr.Dataset({'wh10': wh10,
                             'wh50': wh50})
        print('saving wind ' + year1 + '-' + year2)
        eff_ws.to_netcdf(mer_path+"/eff_ws/merra2_wind_USA_" + year1 + "-" + year2 + ".nc")
        print('calculating alpha ' + year1 + '-' + year2)
        alpha = (xr.ufuncs.log(eff_ws.wh50/eff_ws.wh10)/np.log(50/(10+data.DISPH))).compute()
        print('saving alpha ' + year1 + '-' + year2)
        xr.Dataset({'alpha': alpha}).to_netcdf(mer_path+"/eff_ws/merra2_alpha_USA_" + year1 + "-" + year2 + ".nc")
