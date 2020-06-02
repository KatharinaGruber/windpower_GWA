import xarray as xr
from paths_nz import nz_path

GWA2 = xr.open_rasterio(nz_path[:-2]+'/gwa2_250_ws_DEFLATE.tif')

lat1 = -47.5
lat2 = -34
lon1 = 166
lon2 = 179

GWA2_50 = GWA2.sel(band=1,x=slice(lon1,lon2),y=slice(lat2,lat1))
GWA2_100 = GWA2.sel(band=2,x=slice(lon1,lon2),y=slice(lat2,lat1))

GWA2_50.to_netcdf(nz_path + '/GWA/GWA2_NZ50m.nc')
GWA2_100.to_netcdf(nz_path + '/GWA/GWA2_NZ100m.nc')