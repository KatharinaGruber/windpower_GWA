#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os
import pandas as pd
import xarray as xr
from scipy.stats import pearsonr

from paths_bra import *

from dask.diagnostics import ProgressBar
ProgressBar().register()

# MERRA-2 and ERA5 only unique interpolated locations
print('prepare turbine location data')
# open turbine files
wt_mer = pd.read_csv(bra_path + '/turbine_data_mer.csv', index_col=0)
wt_era = pd.read_csv(bra_path + '/turbine_data_era.csv', index_col=0)

# open wind files
wind_mer = xr.open_mfdataset(mer_path + "/eff_ws/merra2_wind_BRA_*.nc", chunks = {'time': 38})
alpha_mer = xr.open_mfdataset(mer_path + "/eff_ws/merra2_alpha_BRA_*.nc", chunks = {'time': 38})
wind_era = xr.open_mfdataset(era_path + "/eff_ws/era5_wind_BRA_*.nc", chunks = {'time': 38})
alpha_era = xr.open_mfdataset(era_path + "/eff_ws/era5_alpha_BRA_*.nc", chunks = {'time': 38})

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
ip_mer = in_seq_mer.interp(coords={"lon":xr.DataArray(wt_mer.lon,dims='location'),
                                   "lat":xr.DataArray(wt_mer.lat,dims='location')},method="nearest").to_dataframe()
ip_era = in_seq_era.interp(coords={"lon":xr.DataArray(wt_era.lon,dims='location'),
                                   "lat":xr.DataArray(wt_era.lat,dims='location')},method="nearest").to_dataframe()

# find unique locations
uniques_mer = ip_mer.groupby(ip_mer.x).min()
uniques_era = ip_era.groupby(ip_era.x).min()

# add ids to unique locations
uniques_era['cor_id'] = range(len(uniques_era.index))
uniques_mer['cor_id'] = range(len(uniques_mer.index))

# add ids to wind turbine data
wt_mer['cor_id'] = ip_mer.x.map(uniques_mer.cor_id)
wt_era['cor_id'] = ip_era.x.map(uniques_era.cor_id)

ANL = pd.read_csv(bra_path + '/turbine_data.csv', index_col = 0)
lbl = pd.read_csv(bra_path+ '/labels_turbine_data_gwa3.csv',index_col=0)

# Usinas
# some locations have more than one park, get shares of parks
sharesMER = ANL.cap.groupby([lbl.lbl_mer.values,ANL.name.values]).sum()/ANL.cap.groupby([lbl.lbl_mer.values,ANL.name.values]).sum().index.get_level_values(0).map(ANL.cap.groupby(lbl.lbl_mer.values).sum())
sharesERA = ANL.cap.groupby([lbl.lbl_era.values,ANL.name.values]).sum()/ANL.cap.groupby([lbl.lbl_era.values,ANL.name.values]).sum().index.get_level_values(0).map(ANL.cap.groupby(lbl.lbl_era.values).sum())
sharesMERg = ANL.cap.groupby([lbl.lbl_mer_gwa.values,ANL.name.values]).sum()/ANL.cap.groupby([lbl.lbl_mer_gwa.values,ANL.name.values]).sum().index.get_level_values(0).map(ANL.cap.groupby(lbl.lbl_mer_gwa.values).sum())
sharesERAg = ANL.cap.groupby([lbl.lbl_era_gwa.values,ANL.name.values]).sum()/ANL.cap.groupby([lbl.lbl_era_gwa.values,ANL.name.values]).sum().index.get_level_values(0).map(ANL.cap.groupby(lbl.lbl_era_gwa.values).sum())
# add ids to shares
sharesMER = pd.DataFrame({'share':sharesMER,
                          'cor_id':sharesMER.index.get_level_values(0).map(pd.Series(wt_mer.cor_id.values,index=lbl.lbl_mer.unique()))})
sharesERA = pd.DataFrame({'share':sharesERA,
                          'cor_id':sharesERA.index.get_level_values(0).map(pd.Series(wt_era.cor_id.values,index=lbl.lbl_era.unique()))})
# group ids by park
cidMER = sharesMER.groupby(sharesMER.index.get_level_values(1)).cor_id.unique()
cidERA = sharesERA.groupby(sharesERA.index.get_level_values(1)).cor_id.unique()
# get number of grid cells
ngc_USI_mer = cidMER.apply(len)
ngc_USI_era = cidERA.apply(len)

# States
# load matching parks
mpH = pd.read_pickle(bra_path + '/matches2Hlc.pkl')
mpH100 = mpH[mpH.score==100].drop('score',axis=1)
ngc_EST_mer = ngc_USI_mer[mpH100.ANL_name].groupby(mpH100.state.values).sum()
ngc_EST_era = ngc_USI_era[mpH100.ANL_name].groupby(mpH100.state.values).sum()

# Brazil
ngc_BRA_mer = pd.Series(ngc_EST_mer.sum(), index=['BRA'])
ngc_BRA_era = pd.Series(ngc_EST_era.sum(), index=['BRA'])

# merge sizes
print('merge all sizes')
# merge sizes per scale
ngc_USI = pd.concat([pd.DataFrame({'scale':'park',
                                   'region':ngc_USI_mer.index,
                                   'dataset':'MERRA2',
                                   'cor':ngc_USI_mer.values}),
                     pd.DataFrame({'scale':'park',
                                   'region':ngc_USI_era.index,
                                   'dataset':'ERA5',
                                   'cor':ngc_USI_era.values})])

ngc_EST = pd.concat([pd.DataFrame({'scale':'state',
                                   'region':ngc_EST_mer.index,
                                   'dataset':'MERRA2',
                                   'cor':ngc_EST_mer.values}),
                     pd.DataFrame({'scale':'state',
                                   'region':ngc_EST_era.index,
                                   'dataset':'ERA5',
                                   'cor':ngc_EST_era.values})])

ngc_BRA = pd.concat([pd.DataFrame({'scale':'country',
                                   'region':ngc_BRA_mer.index,
                                   'dataset':'MERRA2',
                                   'cor':ngc_BRA_mer.values}),
                     pd.DataFrame({'scale':'country',
                                   'region':ngc_BRA_era.index,
                                   'dataset':'ERA5',
                                   'cor':ngc_BRA_era.values})])
# merge all scales
ngc = pd.concat([ngc_USI,ngc_EST,ngc_BRA]*3, axis=0)
# add temp column (sizes are the same for all)
ngc['temp'] = np.repeat(['m','d','h'],len(ngc)/3) 

# save system sizes calcualted by number of grid points
print('save system sizes')
ngc.to_csv(results_path + '/number_grid_points.csv')