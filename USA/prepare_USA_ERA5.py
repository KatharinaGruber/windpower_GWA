# Script for preparing ERA5 reanalysis data:
# 1. Calculate effective wind speeds from u and v wind speeds in two heights (10 and 100m)
# 2. Calculate alpha friction coefficient

from paths_usa import era_path
import glob
import numpy as np
import xarray as xr

files = glob.glob(era_path+'/*.nc')
files.sort()

out_files = glob.glob(era_path + '/eff_ws/*')

for year in range(0,31,2):
	i1 = year * 12
	i2 = i1 + 24
	wfile = era_path+'/eff_ws/era5_wind_USA_' + str(1987+year) + '-' + str(1987+year+1) + '.nc'
	afile = era_path+'/eff_ws/era5_alpha_USA_' + str(1987+year) + '-' + str(1987+year+1) + '.nc'
	if wfile not in out_files:
		print('calculating wind ' + str(1987+year) + '-' + str(1987+year+1))
		data = xr.open_mfdataset(files[i1:i2], chunks = {'time': 38})
		wh10 = ((data.u10**2+data.v10**2)**0.5).compute()
		wh100 = ((data.u100**2+data.v100**2)**0.5).compute()
		print('saving wind ' + str(1987+year) + '-' + str(1987+year+1))
		eff_ws = xr.Dataset({'wh10': wh10,
							 'wh100': wh100})
						 
		eff_ws.to_netcdf(wfile)
		eff_ws.close()
		del(eff_ws)
	if afile not in out_files:
		print('calculating alpha ' + str(1987+year) + '-' + str(1987+year+1))
		eff_ws = xr.open_dataset(wfile)
		alpha = (xr.ufuncs.log(eff_ws.wh100/eff_ws.wh10)/np.log(100/10)).compute()
		print('saving alpha ' + str(1987+year) + '-' + str(1987+year+1))
		xr.Dataset({'alpha': alpha}).to_netcdf(afile)
		del(alpha)


# 2019
wfile = era_path+'/eff_ws/era5_wind_USA_2019.nc'
afile = era_path+'/eff_ws/era5_alpha_USA_2019.nc'
if wfile not in out_files:
    print('calculating wind 2019')
    data = xr.open_mfdataset(files[-12:], chunks = {'time': 38})
    wh10 = ((data.u10**2+data.v10**2)**0.5).compute()
    wh100 = ((data.u100**2+data.v100**2)**0.5).compute()
    print('saving wind 2019')
    eff_ws = xr.Dataset({'wh10': wh10,
                         'wh100': wh100})
                     
    eff_ws.to_netcdf(wfile)
    eff_ws.close()
    del(eff_ws)
if afile not in out_files:
    print('calculating alpha 2019')
    eff_ws = xr.open_dataset(wfile)
    alpha = (xr.ufuncs.log(eff_ws.wh100/eff_ws.wh10)/np.log(100/10)).compute()
    print('saving alpha 2019')
    xr.Dataset({'alpha': alpha}).to_netcdf(afile)
    del(alpha)