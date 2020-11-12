# Script for preparing ERA5 reanalysis data:
# 1. Calculate effective wind speeds from u and v wind speeds in two heights (10 and 100m)
# 2. Calculate alpha friction coefficient

from paths_nz import era_path
import glob
import numpy as np
import os
import xarray as xr

files = glob.glob(era_path+'/*.nc')
files.sort()

if not os.path.isdir(era_path + '/eff_ws'):
    os.mkdir(era_path + '/eff_ws')

out_files = glob.glob(era_path + '/eff_ws/*')

# 1987 - 2019
i1 = 1987
i2 = 2019
wfile = era_path+'/eff_ws/era5_wind_NZ_' + str(i1) + '-' + str(i2) + '.nc'
afile = era_path+'/eff_ws/era5_alpha_NZ_' + str(i1) + '-' + str(i2) + '.nc'
if wfile not in out_files:
    print('calculating wind ' + str(i1) + '-' + str(i2))
    data = xr.open_mfdataset(files, chunks = {'time': 46})
    wh10 = ((data.u10**2+data.v10**2)**0.5).compute()
    wh100 = ((data.u100**2+data.v100**2)**0.5).compute()
    print('saving wind ' + str(i1) + '-' + str(i2))
    eff_ws = xr.Dataset({'wh10': wh10,
                         'wh100': wh100})
                     
    eff_ws.to_netcdf(wfile)
    eff_ws.close()
    del(eff_ws)
if afile not in out_files:
    print('calculating alpha ' + str(i1) + '-' + str(i2))
    eff_ws = xr.open_dataset(wfile)
    alpha = (xr.ufuncs.log(eff_ws.wh100/eff_ws.wh10)/np.log(100/10)).compute()
    print('saving alpha ' + str(i1) + '-' + str(i2))
    xr.Dataset({'alpha': alpha}).to_netcdf(afile)
    del(alpha)

'''
# split into two periods of 7 years each
for year in [6,13]:
	i1 = (year - 6)*12
	i2 = i1 + 7*12
	wfile = era_path+'/eff_ws/era5_wind_BRA_' + str(2000+year) + '-' + str(2000+year+6) + '.nc'
	afile = era_path+'/eff_ws/era5_alpha_BRA_' + str(2000+year) + '-' + str(2000+year+6) + '.nc'
	if wfile not in out_files:
		print('calculating wind ' + str(2000+year) + '-' + str(2000+year+6))
		data = xr.open_mfdataset(files[i1:i2], chunks = {'time': 60})
		wh10 = ((data.u10**2+data.v10**2)**0.5).compute()
		wh100 = ((data.u100**2+data.v100**2)**0.5).compute()
		print('saving wind ' + str(2000+year) + '-' + str(2000+year+6))
		eff_ws = xr.Dataset({'wh10': wh10,
							 'wh100': wh100})
						 
		eff_ws.to_netcdf(wfile)
		eff_ws.close()
		del(eff_ws)
	if afile not in out_files:
		print('calculating alpha ' + str(2000+year) + '-' + str(2000+year+6))
		eff_ws = xr.open_dataset(wfile)
		alpha = (xr.ufuncs.log(eff_ws.wh100/eff_ws.wh10)/np.log(100/10)).compute()
		print('saving alpha ' + str(2000+year) + '-' + str(2000+year+6))
		xr.Dataset({'alpha': alpha}).to_netcdf(afile)
		del(alpha)
'''