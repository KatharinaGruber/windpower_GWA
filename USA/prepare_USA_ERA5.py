# Script for preparing ERA5 reanalysis data:
# 1. Calculate effective wind speeds from u and v wind speeds in two heights (10 and 100m)
# 2. Calculate alpha friction coefficient

from paths_usa import era_path
import glob

files = glob.glob(era_path+'/*.nc')

out_files = glob.glob(era_path + '/eff_ws/*')

# first period: 2000-2002
wfile = era_path+"/eff_ws/era5_wind_USA_2000-2002.nc"
afile = era_path+"/eff_ws/era5_alpha_USA_2000-2002.nc"
if wfile not in out_files:
	print('calculating wind 2000-2002')
	data = xr.open_mfdataset(files[:25], chunks = {'time': 38})
	wh10 = ((data.u10**2+data.v10**2)**0.5).compute()
	wh100 = ((data.u100**2+data.v100**2)**0.5).compute()

	eff_ws = xr.Dataset({'wh10': wh10,
						 'wh100': wh100})
	print('saving wind 2000-2002')			 
	eff_ws.to_netcdf(wfile)
if afile not in out_files:
	print('calculating alpha 2000-2002')
	eff_ws = xr.open_dataset(wfile)
	alpha = (xr.ufuncs.log(eff_ws.wh100/eff_ws.wh100)/np.log(100/10)).compute()
	print('saving alpha 2000-2002')
	xr.Dataset({'alpha': alpha}).to_netcdf(afile)


# other periods
for year in [3,5,7,9,11,13,15,17]:
	i1 = (year - 1)*12 + 1
	i2 = i1 + 24
	wfile = era_path+'/eff_ws/era5_wind_USA_' + str(2000+year) + '-' + str(2000+year+1) + '.nc'
	afile = era_path+'/eff_ws/era5_alpha_USA_' + str(2000+year) + '-' + str(2000+year+1) + '.nc'
	if wfile not in out_files:
		print('calculating wind ' + str(2000+year) + '-' + str(2000+year+1))
		data = xr.open_mfdataset(files[i1:i2], chunks = {'time': 38})
		wh10 = ((data.u10**2+data.v10**2)**0.5).compute()
		wh100 = ((data.u100**2+data.v100**2)**0.5).compute()
		print('saving wind ' + str(2000+year) + '-' + str(2000+year+1))
		eff_ws = xr.Dataset({'wh10': wh10,
							 'wh100': wh100})
						 
		eff_ws.to_netcdf(wfile)
	if afile not in out_files:
		print('calculating alpha ' + str(2000+year) + '-' + str(2000+year+1))
		eff_ws = xr.open_dataset(wfile)
		alpha = (xr.ufuncs.log(eff_ws.wh100/eff_ws.wh10)/np.log(100/10)).compute()
		print('saving alpha ' + str(2000+year) + '-' + str(2000+year+1))
		xr.Dataset({'alpha': alpha}).to_netcdf(afile)
