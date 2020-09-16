import glob
import os
import xarray as xr
from paths_bra import bra_path


if bra_path + '/GWA/GWA2_BRA50m.nc' not in glob.glob(bra_path + '/GWA/*.nc'):
    if not os.path.isdir(bra_path + '/GWA'):
        os.mkdir(bra_path + '/GWA')
    GWA2 = xr.open_rasterio(bra_path[:-3]+'/gwa2_250_ws_DEFLATE.tif')

    lat1 = -36
    lat2 = 5.5
    lon1 = -74.1
    lon2 = -33
    
    GWA2_50 = GWA2.sel(band=1,x=slice(lon1,lon2),y=slice(lat2,lat1))
    GWA2_100 = GWA2.sel(band=2,x=slice(lon1,lon2),y=slice(lat2,lat1))
    GWA2_200 = GWA2.sel(band=3,x=slice(lon1,lon2),y=slice(lat2,lat1))

    GWA2_50.to_netcdf(bra_path + '/GWA/GWA2_BRA50m.nc')
    GWA2_100.to_netcdf(bra_path + '/GWA/GWA2_BRA100m.nc')
    GWA2_200.to_netcdf(bra_path + '/GWA/GWA2_BRA200m.nc')