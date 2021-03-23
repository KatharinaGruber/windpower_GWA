import glob
import numpy as np
import xarray as xr

from paths import *

startGWA2 = '1987'
startGWA3 = '2008'
endGWA2 = '2016'
endGWA3 = '2017'


# South Africa
print('South Africa')
# ERA5
windE_ZAF = xr.open_dataset(era_path + '/ZAF/eff_ws/era5_wind_ZAF_1987-2019.nc')
# GWA2
if results_path + '/ZAF/cf_ERA_GWA2.nc' not in glob.glob(results_path + '/ZAF/*.nc'):
    print('South Africa ERA5 GWA2')
    GWA2E_ZAF = xr.open_dataarray(data_path + '/ZAF/GWA/GWA2_ZAF100m.nc').drop('band')
    mwindEg2_ZAF = windE_ZAF.wh100.sel(time=slice(startGWA2,endGWA2)
                                      ).mean('time').interp(coords={'longitude':GWA2E_ZAF.x.values,
                                                                    'latitude':GWA2E_ZAF.y.values},
                                                            method = 'nearest')
    cfEg2_ZAF = GWA2E_ZAF.where(GWA2E_ZAF>=0,np.nan).values/mwindEg2_ZAF
    cfEg2_ZAF.name = 'cf'
    cfEg2_ZAF.to_netcdf(results_path + '/ZAF/cf_ERA_GWA2.nc')
    del(mwindEg2_ZAF,GWA2E_ZAF,cfEg2_ZAF)                                                                
    print('South Africa ERA5 GWA2 done')
# GWA3
if results_path + '/ZAF/cf_ERA_GWA3.nc' not in glob.glob(results_path + '/ZAF/*.nc'):
    print('South Africa ERA5 GWA3')
    GWA3E_ZAF = xr.open_rasterio(data_path + '/ZAF/GWA/GWA3_ZAF100m.tif').sel(band=1).drop('band').dropna(dim='x',how='all').dropna(dim='y',how='all')
    mwindEg3_ZAF = windE_ZAF.wh100.sel(time=slice(startGWA3,endGWA3)
                                      ).mean('time').interp(coords={'longitude':GWA3E_ZAF.x.values,
                                                                    'latitude':GWA3E_ZAF.y.values},
                                                            method = 'nearest')
    cfEg3_ZAF = GWA3E_ZAF.where(GWA3E_ZAF>=0,np.nan).values/mwindEg3_ZAF
    cfEg3_ZAF.name = 'cf'
    cfEg3_ZAF.to_netcdf(results_path + '/ZAF/cf_ERA_GWA3.nc')
    del(mwindEg3_ZAF,GWA3E_ZAF,cfEg3_ZAF)
    print('South Africa ERA5 GWA3 done')
del(windE_ZAF)
# MERRA-2
windM_ZAF = xr.open_dataset(mer_path + '/ZAF/eff_ws/merra2_wind_ZAF_1987-2019.nc')
# GWA2
if results_path + '/ZAF/cf_MERRA_GWA2.nc' not in glob.glob(results_path + '/ZAF/*.nc'):
    print('South Africa MERRA-2 GWA2')
    GWA2M_ZAF = xr.open_dataarray(data_path + '/ZAF/GWA/GWA2_ZAF50m.nc').drop('band')
    mwindMg2_ZAF = windM_ZAF.wh50.sel(time=slice(startGWA2,endGWA2)
                                      ).mean('time').interp(coords={'lon':GWA2M_ZAF.x.values,
                                                                    'lat':GWA2M_ZAF.y.values},
                                                            method = 'nearest')
    cfMg2_ZAF = GWA2M_ZAF.where(GWA2M_ZAF>=0,np.nan).values/mwindMg2_ZAF
    cfMg2_ZAF.name = 'cf'
    cfMg2_ZAF.to_netcdf(results_path + '/ZAF/cf_MERRA_GWA2.nc')
    del(mwindMg2_ZAF,GWA2M_ZAF,cfMg2_ZAF)
    print('South Africa MERRA-2 GWA2 done')
# GWA3
if results_path + '/ZAF/cf_MERRA_GWA3.nc' not in glob.glob(results_path + '/ZAF/*.nc'):
    print('South Africa MERRA-2 GWA3')
    GWA3M_ZAF = xr.open_rasterio(data_path + '/ZAF/GWA/GWA3_ZAF50m.tif').sel(band=1).drop('band').dropna(dim='x',how='all').dropna(dim='y',how='all')
    mwindMg3_ZAF = windM_ZAF.wh50.sel(time=slice(startGWA3,endGWA3)
                                      ).mean('time').interp(coords={'lon':GWA3M_ZAF.x.values,
                                                                    'lat':GWA3M_ZAF.y.values},
                                                            method = 'nearest')
    cfMg3_ZAF = GWA3M_ZAF.where(GWA3M_ZAF>=0,np.nan).values/mwindMg3_ZAF
    cfMg3_ZAF.name = 'cf'
    cfMg3_ZAF.to_netcdf(results_path + '/ZAF/cf_MERRA_GWA3.nc')
    del(mwindMg3_ZAF,GWA3M_ZAF,cfMg3_ZAF)
    print('South Africa MERRA-2 GWA3 done')
del(windM_ZAF)

# New Zealand
print('New Zealand')
# ERA5
windE_NZ = xr.open_dataset(era_path + '/NZ/eff_ws/era5_wind_NZ_1987-2019.nc')
# GWA2
if results_path + '/NZ/cf_ERA_GWA2.nc' not in glob.glob(results_path + '/NZ/*.nc'):
    print('New Zealand ERA5 GWA2')
    GWA2E_NZ = xr.open_dataarray(data_path + '/NZ/GWA/GWA2_NZ100m.nc').drop('band')
    mwindEg2_NZ = windE_NZ.wh100.sel(time=slice(startGWA2,endGWA2)
                                      ).mean('time').interp(coords={'longitude':GWA2E_NZ.x.values,
                                                                    'latitude':GWA2E_NZ.y.values},
                                                            method = 'nearest')
    cfEg2_NZ = GWA2E_NZ.where(GWA2E_NZ>=0,np.nan).values/mwindEg2_NZ
    cfEg2_NZ.name = 'cf'
    cfEg2_NZ.to_netcdf(results_path + '/NZ/cf_ERA_GWA2.nc')
    del(mwindEg2_NZ,GWA2E_NZ,cfEg2_NZ)
    print('New Zealand ERA5 GWA2 done')
# GWA3
if results_path + '/NZ/cf_ERA_GWA3.nc' not in glob.glob(results_path + '/NZ/*.nc'):
    print('New Zealand ERA5 GWA3')
    GWA3E_NZ = xr.open_rasterio(data_path + '/NZ/GWA/GWA3_NZ100m.tif').sel(band=1).drop('band').dropna(dim='x',how='all').dropna(dim='y',how='all')
    mwindEg3_NZ = windE_NZ.wh100.sel(time=slice(startGWA3,endGWA3)
                                      ).mean('time').interp(coords={'longitude':GWA3E_NZ.x.values,
                                                                    'latitude':GWA3E_NZ.y.values},
                                                            method = 'nearest')
    cfEg3_NZ = GWA3E_NZ.where(GWA3E_NZ>=0,np.nan).values/mwindEg3_NZ
    cfEg3_NZ.name = 'cf'
    cfEg3_NZ.to_netcdf(results_path + '/NZ/cf_ERA_GWA3.nc')
    del(mwindEg3_NZ,GWA3E_NZ,cfEg3_NZ)
    print('New Zealand ERA5 GWA3 done')
del(windE_NZ)
# MERRA-2
windM_NZ = xr.open_dataset(mer_path + '/NZ/eff_ws/merra2_wind_NZ_1987-2019.nc')
# GWA2
if results_path + '/NZ/cf_MERRA_GWA2.nc' not in glob.glob(results_path + '/NZ/*.nc'):
    print('New Zealand MERRA-2 GWA2')
    GWA2M_NZ = xr.open_dataarray(data_path + '/NZ/GWA/GWA2_NZ50m.nc').drop('band')
    mwindMg2_NZ = windM_NZ.wh50.sel(time=slice(startGWA2,endGWA2)
                                    ).mean('time').interp(coords={'lon':GWA2M_NZ.x.values,
                                                                    'lat':GWA2M_NZ.y.values},
                                                            method = 'nearest')
    cfMg2_NZ = GWA2M_NZ.where(GWA2M_NZ>=0,np.nan).values/mwindMg2_NZ
    cfMg2_NZ.name = 'cf'
    cfMg2_NZ.to_netcdf(results_path + '/NZ/cf_MERRA_GWA2.nc')
    del(mwindMg2_NZ,GWA2M_NZ,cfMg2_NZ)
    print('New Zealand MERRA-2 GWA2 done')
# GWA3
if results_path + '/NZ/cf_MERRA_GWA3.nc' not in glob.glob(results_path + '/NZ/*.nc'):
    print('New Zealand MERRA-2 GWA3')
    GWA3M_NZ = xr.open_rasterio(data_path + '/NZ/GWA/GWA3_NZ50m.tif').sel(band=1).drop('band').dropna(dim='x',how='all').dropna(dim='y',how='all')
    mwindMg3_NZ = windM_NZ.wh50.sel(time=slice(startGWA3,endGWA3)
                                      ).mean('time').interp(coords={'lon':GWA3M_NZ.x.values,
                                                                    'lat':GWA3M_NZ.y.values},
                                                            method = 'nearest')
    cfMg3_NZ = GWA3M_NZ.where(GWA3M_NZ>=0,np.nan).values/mwindMg3_NZ
    cfMg3_NZ.name = 'cf'           
    cfMg3_NZ.to_netcdf(results_path + '/NZ/cf_MERRA_GWA3.nc')
    del(mwindMg3_NZ,GWA3M_NZ,cfMg3_NZ)
    print('New Zealand MERRA-2 GWA3 done')
del(windM_NZ)                                                                
                                                                
# Brazil
print('Brazil')
# ERA5
windE_BRA = xr.open_mfdataset(era_path + '/BRA/eff_ws/era5_wind_BRA_*.nc')                                                             
# GWA2
if results_path + '/BRA/cf_ERA_GWA2.nc' not in glob.glob(results_path + '/BRA/*.nc'):
    print('Brazil ERA5 GWA2')
    GWA2E_BRA = xr.open_dataarray(data_path + '/BRA/GWA/GWA2_BRA100m.nc').drop('band')
    mwindEg2_BRA = windE_BRA.wh100.sel(time=slice(startGWA2,endGWA2)
                                      ).mean('time').interp(coords={'longitude':GWA2E_BRA.x.values,
                                                                    'latitude':GWA2E_BRA.y.values},
                                                            method = 'nearest')
    cfEg2_BRA = GWA2E_BRA.where(GWA2E_BRA>=0,np.nan).values/mwindEg2_BRA
    cfEg2_BRA.name = 'cf'
    cfEg2_BRA.to_netcdf(results_path + '/BRA/cf_ERA_GWA2.nc')
    del(mwindEg2_BRA,GWA2E_BRA,cfEg2_BRA)
    print('Brazil ERA5 GWA2 done')
# GWA3
if results_path + '/BRA/cf_ERA_GWA3.nc' not in glob.glob(results_path + '/BRA/*.nc'):
    print('Brazil ERA5 GWA3')
    GWA3E_BRA = xr.open_rasterio(data_path + '/BRA/GWA/GWA3_BRA100m.tif').sel(band=1).drop('band').dropna(dim='x',how='all').dropna(dim='y',how='all')
    mwindEg3_BRA = windE_BRA.wh100.sel(time=slice(startGWA3,endGWA3)
                                      ).mean('time').interp(coords={'longitude':GWA3E_BRA.x.values,
                                                                    'latitude':GWA3E_BRA.y.values},
                                                            method = 'nearest')
    cfEg3_BRA = GWA3E_BRA.where(GWA3E_BRA>=0,np.nan).values/mwindEg3_BRA
    cfEg3_BRA.name = 'cf'
    cfEg3_BRA.to_netcdf(results_path + '/BRA/cf_ERA_GWA3.nc')
    del(mwindEg3_BRA,GWA3E_BRA,cfEg3_BRA)
    print('Brazil ERA5 GWA3 done')
del(windE_BRA)
# MERRA-2
windM_BRA = xr.open_mfdataset(mer_path + '/BRA/eff_ws/merra2_wind_BRA_*.nc')
# GWA2
if results_path + '/BRA/cf_MERRA_GWA2.nc' not in glob.glob(results_path + '/BRA/*.nc'):
    print('Brazil MERRA-2 GWA2')
    GWA2M_BRA = xr.open_dataarray(data_path + '/BRA/GWA/GWA2_BRA50m.nc').drop('band')
    mwindMg2_BRA = windM_BRA.wh50.sel(time=slice(startGWA2,endGWA2)
                                      ).mean('time').interp(coords={'lon':GWA2M_BRA.x.values,
                                                                    'lat':GWA2M_BRA.y.values},
                                                            method = 'nearest')
    cfMg2_BRA = GWA2M_BRA.where(GWA2M_BRA>=0,np.nan).values/mwindMg2_BRA
    cfMg2_BRA.name = 'cf'
    cfMg2_BRA.to_netcdf(results_path + '/BRA/cf_MERRA_GWA2.nc')
    del(mwindMg2_BRA,GWA2M_BRA,cfMg2_BRA)
    print('Brazil MERRA-2 GWA2 done')
# GWA3
if results_path + '/BRA/cf_MERRA_GWA3.nc' not in glob.glob(results_path + '/BRA/*.nc'):
    print('Brazil MERRA-2 GWA3')
    GWA3M_BRA = xr.open_rasterio(data_path + '/BRA/GWA/GWA3_BRA50m.tif').sel(band=1).drop('band').dropna(dim='x',how='all').dropna(dim='y',how='all')
    mwindMg3_BRA = windM_BRA.wh50.sel(time=slice(startGWA3,endGWA3)
                                      ).mean('time').interp(coords={'lon':GWA3M_BRA.x.values,
                                                                    'lat':GWA3M_BRA.y.values},
                                                            method = 'nearest')                                                           
    cfMg3_BRA = GWA3M_BRA.where(GWA3M_BRA>=0,np.nan).values/mwindMg3_BRA
    cfMg3_BRA.name = 'cf'
    cfMg3_BRA.to_netcdf(results_path + '/BRA/cf_MERRA_GWA3.nc')
    del(mwindMg3_BRA,GWA3M_BRA,cfMg3_BRA)
    print('Brazil MERRA-2 GWA3 done')
del(windM_BRA)                                 

# USA
print('USA')
# ERA5
windE_USA = xr.open_mfdataset(era_path + '/USA/eff_ws/era5_wind_USA_*.nc')
# GWA2
if results_path + '/USA/cf_ERA_GWA2.nc' not in glob.glob(results_path + '/USA/*.nc'):
    print('USA ERA5 GWA2')
    GWA2E_USA = xr.open_rasterio(data_path + '/USA/GWA/GWA_USA100m.tif').sel(band=1).drop('band').dropna(dim='x',how='all').dropna(dim='y',how='all')
    latmin = max(min(windE_USA.latitude.values),min(GWA2E_USA.y.values))
    latmax = min(max(windE_USA.latitude.values),max(GWA2E_USA.y.values))
    lonmin = max(min(windE_USA.longitude.values),min(GWA2E_USA.x.values))
    lonmax = min(max(windE_USA.longitude.values),max(GWA2E_USA.x.values))
    GWA2E_USAcut = GWA2E_USA.sel(x=slice(lonmin,lonmax),y=slice(latmax,latmin))    
    mwindEg2_USA = windE_USA.wh100.sel(time=slice(startGWA2,endGWA2)
                                      ).mean('time').interp(coords={'longitude':GWA2E_USAcut.x.values,
                                                                    'latitude':GWA2E_USAcut.y.values},
                                                            method = 'nearest')
    cfEg2_USA = GWA2E_USAcut.where(GWA2E_USAcut>=0,np.nan).values/mwindEg2_USA
    cfEg2_USA.name = 'cf'
    cfEg2_USA.to_netcdf(results_path + '/USA/cf_ERA_GWA2.nc')
    del(mwindEg2_USA,GWA2E_USA,GWA2E_USAcut,cfEg2_USA)
    print('USA ERA5 GWA2 done')
# GWA3
if results_path + '/USA/cf_ERA_GWA3.nc' not in glob.glob(results_path + '/USA/*.nc'):
    print('USA ERA5 GWA3')
    GWA3E_USA = xr.open_rasterio(data_path + '/USA/GWA/GWA3_USA100m.tif').sel(band=1).drop('band').dropna(dim='x',how='all').dropna(dim='y',how='all')
    latmin = max(min(windE_USA.latitude.values),min(GWA3E_USA.y.values))
    latmax = min(max(windE_USA.latitude.values),max(GWA3E_USA.y.values))
    lonmin = max(min(windE_USA.longitude.values),min(GWA3E_USA.x.values))
    lonmax = min(max(windE_USA.longitude.values),max(GWA3E_USA.x.values))
    GWA3E_USAcut = GWA3E_USA.sel(x=slice(lonmin,lonmax),y=slice(latmax,latmin))    
    mwindEg3_USA = windE_USA.wh100.sel(time=slice(startGWA3,endGWA3)
                                      ).mean('time').interp(coords={'longitude':GWA3E_USAcut.x.values,
                                                                    'latitude':GWA3E_USAcut.y.values},
                                                            method = 'nearest')
    cfEg3_USA = GWA3E_USAcut.where(GWA3E_USAcut>=0,np.nan).values/mwindEg3_USA
    cfEg3_USA.name = 'cf'
    cfEg3_USA.to_netcdf(results_path + '/USA/cf_ERA_GWA3.nc')
    del(mwindEg3_USA,GWA3E_USA,GWA3E_USAcut,cfEg3_USA)
    print('USA ERA5 GWA3 done')
del(windE_USA)                                                                
# MERRA-2
windM_USA = xr.open_mfdataset(mer_path + '/USA/eff_ws/merra2_wind_USA_*.nc')
# GWA2
if results_path + '/USA/cf_MERRA_GWA2.nc' not in glob.glob(results_path + '/USA/*.nc'):
    print('USA MERRA-2 GWA2')
    GWA2M_USA = xr.open_rasterio(data_path + '/USA/GWA/GWA_USA50m.tif').sel(band=1).drop('band').dropna(dim='x',how='all').dropna(dim='y',how='all')
    latmin = max(min(windM_USA.lat.values),min(GWA2M_USA.y.values))
    latmax = min(max(windM_USA.lat.values),max(GWA2M_USA.y.values))
    lonmin = max(min(windM_USA.lon.values),min(GWA2M_USA.x.values))
    lonmax = min(max(windM_USA.lon.values),max(GWA2M_USA.x.values))
    GWA2M_USAcut = GWA2M_USA.sel(x=slice(lonmin,lonmax),y=slice(latmax,latmin))
    mwindMg2_USA = windM_USA.wh50.sel(time=slice(startGWA2,endGWA2)
                                      ).mean('time').interp(coords={'lon':GWA2M_USAcut.x.values,
                                                                    'lat':GWA2M_USAcut.y.values},
                                                            method = 'nearest')
    cfMg2_USA = GWA2M_USAcut.where(GWA2M_USAcut>=0,np.nan).values/mwindMg2_USA
    cfMg2_USA.name = 'cf'
    cfMg2_USA.to_netcdf(results_path + '/USA/cf_MERRA_GWA2.nc')
    del(mwindMg2_USA,GWA2M_USA,GWA2M_USAcut,cfMg2_USA)
    print('USA MERRA-2 GWA2 done')
# GWA3
if results_path + '/USA/cf_MERRA_GWA3.nc' not in glob.glob(results_path + '/USA/*.nc'):
    print('USA MERRA-2 GWA3')
    GWA3M_USA = xr.open_rasterio(data_path + '/USA/GWA/GWA3_USA50m.tif').sel(band=1).drop('band').dropna(dim='x',how='all').dropna(dim='y',how='all')
    latmin = max(min(windM_USA.lat.values),min(GWA3M_USA.y.values))
    latmax = min(max(windM_USA.lat.values),max(GWA3M_USA.y.values))
    lonmin = max(min(windM_USA.lon.values),min(GWA3M_USA.x.values))
    lonmax = min(max(windM_USA.lon.values),max(GWA3M_USA.x.values))
    GWA3M_USAcut = GWA3M_USA.sel(x=slice(lonmin,lonmax),y=slice(latmax,latmin))
    mwindMg3_USA = windM_USA.wh50.sel(time=slice(startGWA3,endGWA3)
                                      ).mean('time').interp(coords={'lon':GWA3M_USAcut.x.values,
                                                                    'lat':GWA3M_USAcut.y.values},
                                                            method = 'nearest')
    cfMg3_USA = GWA3M_USAcut.where(GWA3M_USAcut>=0,np.nan).values/mwindMg3_USA
    cfMg3_USA.name = 'cf'
    cfMg3_USA.to_netcdf(results_path + '/USA/cf_MERRA_GWA3.nc')
    del(mwindMg3_USA,GWA3M_USA,GWA3M_USAcut,cfMg3_USA)
    print('USA MERRA-2 GWA3 done')
del(windM_USA)
