#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os
import pandas as pd
import xarray as xr

from paths_zaf import *

from dask.diagnostics import ProgressBar
ProgressBar().register()

# load windpark data
windparks = pd.read_csv(zaf_path + "/windparks_ZAF.csv", parse_dates=['commissioning'])

# open wind files
wind_mer = xr.open_dataset(mer_path + "/eff_ws/merra2_wind_ZAF_1987-2019.nc", chunks = {'time': 100})
alpha_mer = xr.open_dataset(mer_path + "/eff_ws/merra2_alpha_ZAF_1987-2019.nc", chunks = {'time': 100})
wind_era = xr.open_dataset(era_path + "/eff_ws/era5_wind_ZAF_1987-2019.nc", chunks = {'time': 100})
alpha_era = xr.open_dataset(era_path + "/eff_ws/era5_alpha_ZAF_1987-2019.nc", chunks = {'time': 100})

# Create dataframe with sequence the size of MERRA-2 grid to find out which turbines interpolate to the same point
in_seq_mer = xr.Dataset({'x':(['lat','lon'],
                              np.array(range(wind_mer.wh50.isel(time=0).values.size)).reshape(wind_mer.wh50.isel(time=0).values.shape))},
                         coords = {'lat':wind_mer.lat.values,
                                   'lon':wind_mer.lon.values})
in_seq_era = xr.Dataset({'x':(['lat','lon'],
                              np.array(range(wind_era.wh100.isel(time=0).values.size)).reshape(wind_era.wh100.isel(time=0).values.shape))},
                         coords = {'lat':wind_era.latitude.values,
                                   'lon':wind_era.longitude.values})

# interpolate to reanalysis grid points
ip_mer = in_seq_mer.interp(coords={"lon":xr.DataArray(windparks.Longitude,dims='location'),
                                   "lat":xr.DataArray(windparks.Latitude,dims='location')},method="nearest").to_dataframe()
ip_era = in_seq_era.interp(coords={"lon":xr.DataArray(windparks.Longitude,dims='location'),
                                   "lat":xr.DataArray(windparks.Latitude,dims='location')},method="nearest").to_dataframe()

# merge sizes
print('merge all sizes')
# merge sizes per scale
ngc_cape = pd.DataFrame({'scale':'state',
                         'dataset':np.repeat(['ERA5','MERRA2'],len(windparks.Area.unique())),
                         'cor':pd.concat([ip_era.groupby(windparks.Area.rename('region')).apply(lambda p: len(p.x.unique())),
                                          ip_era.groupby(windparks.Area.rename('region')).apply(lambda p: len(p.x.unique()))],axis=0)}).reset_index()
ngc_ZAF = pd.DataFrame({'region':'ZAF',
                        'scale':'country',
                        'dataset':['ERA5','MERRA2'],
                        'cor':[ip_era.x.unique().shape[0],ip_mer.x.unique().shape[0]]})
# merge all scales
ngc = pd.concat([ngc_cape,ngc_ZAF]*3, axis=0)
# add temp column (sizes are the same for all)
ngc['temp'] = np.repeat(['m','d','h'],len(ngc)/3)
# save system sizes calcualted by number of grid points
print('save system sizes')
ngc.to_csv(results_path + '/number_grid_points.csv')