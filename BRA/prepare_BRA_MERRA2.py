# Script for preparing MERRA-2 reanalysis data:
# 1. Calculate effective wind speeds from u and v wind speeds in two heights (10 and 50m)
# 2. Calculate alpha friction coefficient

from paths_bra import mer_path
import glob
import numpy as np
import os
import xarray as xr

files = glob.glob(mer_path+'/*.nc')
files.sort()

if not os.path.isdir(mer_path + '/eff_ws'):
    os.mkdir(mer_path + '/eff_ws')

out_files = glob.glob(mer_path + '/eff_ws/*')

# 1987 - 1991
i1 = 0
i2 = 5*12
wfile = mer_path+'/eff_ws/merra2_wind_BRA_1987-1991.nc'
afile = mer_path+'/eff_ws/merra2_alpha_BRA_1987-1991.nc'
if wfile not in out_files:
    print('calculating wind 1987-1991')
    data = xr.open_mfdataset(files[i1:i2], chunks = {'time': 60})
    wh10 = ((data.U10M**2+data.V10M**2)**0.5).compute()
    wh50 = ((data.U50M**2+data.V50M**2)**0.5).compute()
    print('saving wind 1987-1991')
    eff_ws = xr.Dataset({'wh10': wh10,
                         'wh50': wh50})
    eff_ws.to_netcdf(wfile)
    eff_ws.close()
    del(eff_ws)
if afile not in out_files:
    print('calculating alpha 1987-1991')
    eff_ws = xr.open_dataset(wfile)
    alpha = (xr.ufuncs.log(eff_ws.wh50/eff_ws.wh10)/np.log(50/10)).compute()
    print('saving alpha 1987-1991')
    xr.Dataset({'alpha': alpha}).to_netcdf(afile)
    del(alpha)

# 1992 - 2019
# split into four periods of 7 years each
for year in range(0,28,7):
	i1 = year*12 + 5*12
	i2 = i1 + 7*12
	wfile = mer_path+'/eff_ws/merra2_wind_BRA_' + str(1992+year) + '-' + str(1992+year+6) + '.nc'
	afile = mer_path+'/eff_ws/merra2_alpha_BRA_' + str(1992+year) + '-' + str(1992+year+6) + '.nc'
	if wfile not in out_files:
		print('calculating wind ' + str(1992+year) + '-' + str(1992+year+6))
		data = xr.open_mfdataset(files[i1:i2], chunks = {'time': 60})
		wh10 = ((data.U10M**2+data.V10M**2)**0.5).compute()
		wh50 = ((data.U50M**2+data.V50M**2)**0.5).compute()
		print('saving wind ' + str(1992+year) + '-' + str(1992+year+6))
		eff_ws = xr.Dataset({'wh10': wh10,
							 'wh50': wh50})
						 
		eff_ws.to_netcdf(wfile)
		eff_ws.close()
		del(eff_ws)
	if afile not in out_files:
		print('calculating alpha ' + str(1992+year) + '-' + str(1992+year+6))
		eff_ws = xr.open_dataset(wfile)
		alpha = (xr.ufuncs.log(eff_ws.wh50/eff_ws.wh10)/np.log(50/(10+data.DISPH))).compute()
		print('saving alpha ' + str(1992+year) + '-' + str(1992+year+6))
		xr.Dataset({'alpha': alpha}).to_netcdf(afile)
		del(alpha)
