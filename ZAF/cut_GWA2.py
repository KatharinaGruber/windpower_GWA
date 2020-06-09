import glob
import os
import xarray as xr
from paths_zaf import zaf_path


if zaf_path + '/GWA/GWA2_ZAF50m.nc' not in glob.glob(zaf_path + '/*.nc'):
    if not os.path.isdir(zaf_path + '/GWA'):
        os.mkdir(zaf_path + '/GWA')
    GWA2 = xr.open_rasterio(zaf_path[:-3]+'/gwa2_250_ws_DEFLATE.tif')

    lat1 = -35
    lat2 = -22
    lon1 = 16
    lon2 = 33

    GWA2_50 = GWA2.sel(band=1,x=slice(lon1,lon2),y=slice(lat2,lat1))
    GWA2_100 = GWA2.sel(band=2,x=slice(lon1,lon2),y=slice(lat2,lat1))

    GWA2_50.to_netcdf(zaf_path + '/GWA/GWA2_ZAF50m.nc')
    GWA2_100.to_netcdf(zaf_path + '/GWA/GWA2_ZAF100m.nc')