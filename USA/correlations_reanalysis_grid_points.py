#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os
import pandas as pd
import xarray as xr
from scipy.stats import pearsonr

from paths_usa import *

from dask.diagnostics import ProgressBar
ProgressBar().register()


# MERRA-2 and ERA5 only unique interpolated locations
print('prepare turbine location data')
# open turbine files
wt_mer = pd.read_csv(usa_path + '/turbine_data_mer.csv', index_col=0)
wt_era = pd.read_csv(usa_path + '/turbine_data_era.csv', index_col=0)

# open wind files
wind_mer = xr.open_mfdataset(mer_path + "/eff_ws/merra2_wind_USA_*.nc", chunks = {'time': 38})
alpha_mer = xr.open_mfdataset(mer_path + "/eff_ws/merra2_alpha_USA_*.nc", chunks = {'time': 38})
wind_era = xr.open_mfdataset(era_path + "/eff_ws/era5_wind_USA_*.nc", chunks = {'time': 38})
alpha_era = xr.open_mfdataset(era_path + "/eff_ws/era5_alpha_USA_*.nc", chunks = {'time': 38})

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

# add ids to unique correlation locations
uniques_era['cor_id'] = range(len(uniques_era.index))
uniques_mer['cor_id'] = range(len(uniques_mer.index))

# add correlation ids to wind turbine data
wt_mer['cor_id'] = ip_mer.x.map(uniques_mer.cor_id)
wt_era['cor_id'] = ip_era.x.map(uniques_era.cor_id)

# interpolate wind to unique locations and extrapolate to 80 m hubheight, which is mean hubheight in US
windi_mer = wind_mer.interp(coords={"lon":xr.DataArray(uniques_mer.lon,dims='location'),
                                    "lat":xr.DataArray(uniques_mer.lat,dims='location')},method="nearest")
windi_era = wind_era.interp(coords={"longitude":xr.DataArray(uniques_era.lon,dims='location'),
                                    "latitude":xr.DataArray(uniques_era.lat,dims='location')},method="nearest")
                                    
alphai_mer = alpha_mer.interp(coords={"lon":xr.DataArray(uniques_mer.lon,dims='location'),
                                      "lat":xr.DataArray(uniques_mer.lat,dims='location')},method="nearest")
alphai_era = alpha_era.interp(coords={"longitude":xr.DataArray(uniques_era.lon,dims='location'),
                                      "latitude":xr.DataArray(uniques_era.lat,dims='location')},method="nearest")
# calculate wind speeds at 80 m height
windhh_mer = (windi_mer.wh50 * (80/50)**alphai_mer.alpha)
windhh_era = (windi_era.wh100 * (80/100)**alphai_era.alpha)


# calculate cross-correlations
print('calculate cross-correlations')
if not os.path.isfile(results_path + '/corr_mer.nc'):
    c_mer = np.corrcoef(windhh_mer,rowvar=False)
    xr.DataArray(c_mer,
                 dims = ['id1','id2'],
                 coords = {'id1':range(c_mer.shape[0]),
                           'id2':range(c_mer.shape[0])}).to_netcdf(results_path + '/corr_mer.nc')
    del(c_mer)
if not os.path.isfile(results_path + '/corr_era.nc'):
    c_era = np.corrcoef(windhh_era,rowvar=False)
    xr.DataArray(c_era,
                 dims = ['id1','id2'],
                 coords = {'id1':range(c_era.shape[0]),
                           'id2':range(c_era.shape[0])}).to_netcdf(results_path + '/corr_era.nc')
    del(c_era)

# load correlations
c_mer = xr.open_dataarray(results_path + '/corr_mer.nc').values
c_era = xr.open_dataarray(results_path + '/corr_era.nc').values

# BPA
print('correlations BPA')

# get BPA windparks
BPA_parks = pd.read_csv(usa_path + "/BPA_windparks.csv")
# get windturbine locations/names
windturbines = pd.read_csv(usa_path + "/uswtdb_v2_3_20200109.csv",delimiter=',')
# get labels
labels = pd.read_csv(usa_path + '/labels_turbine_data_gwa3.csv')
# get indices of BPA wind parks from wind turbine dataset
pBPA = pd.DataFrame({'p': [park in BPA_parks.name.values for park in windturbines[windturbines.t_state!='GU'].p_name.values]})

# remove turbines that are in other states - error in data?
pBPA[windturbines[windturbines.t_state!='GU'].xlong.values<-125] = False
pBPA[windturbines[windturbines.t_state!='GU'].xlong.values>-115] = False

# get shares of where BPA wind parks are installed
BPA_mer = pBPA.groupby(labels.lbl_mer).mean()
BPA_era = pBPA.groupby(labels.lbl_era).mean()
# get correlation IDs
BPA_mer['cor_id'] = wt_mer.cor_id.values
BPA_era['cor_id'] = wt_era.cor_id.values
# select only windparks in BPA
BPA_mer = BPA_mer[BPA_mer.p>0]
BPA_era = BPA_era[BPA_era.p>0]

# unique mean correlations BPA
mc_BPA_mer = np.array([list(c_mer[BPA_mer.cor_id.unique()[i],
                                 np.delete(BPA_mer.cor_id.unique(),i,0)]) for i in range(len(BPA_mer.cor_id.unique()))]).mean()
mc_BPA_era = np.array([list(c_era[BPA_era.cor_id.unique()[i],
                                 np.delete(BPA_era.cor_id.unique(),i,0)]) for i in range(len(BPA_era.cor_id.unique()))]).mean()

# USA
print('correlations USA')

# unique mean correlations USA
mc_USA_mer = np.array([list(c_mer[wt_mer.cor_id.unique()[i],
                                np.delete(wt_mer.cor_id.unique(),i,0)]) for i in range(len(wt_mer.cor_id.unique()))]).mean()
mc_USA_era = np.array([list(c_era[wt_era.cor_id.unique()[i],
                                np.delete(wt_era.cor_id.unique(),i,0)]) for i in range(len(wt_era.cor_id.unique()))]).mean()


# New England
print('correlations New England')

NE_states = ['CT','NH','ME','MA','RI','VT']
wt_mer_NE = wt_mer[[state in NE_states for state in wt_mer.state]]
wt_era_NE = wt_era[[state in NE_states for state in wt_era.state]]
# unique mean correlations New England
mc_NE_mer = np.array([list(c_mer[wt_mer_NE.cor_id.unique()[i],
                                np.delete(wt_mer_NE.cor_id.unique(),i,0)]) for i in range(len(wt_mer_NE.cor_id.unique()))]).mean()
mc_NE_era = np.array([list(c_era[wt_era_NE.cor_id.unique()[i],
                                np.delete(wt_era_NE.cor_id.unique(),i,0)]) for i in range(len(wt_era_NE.cor_id.unique()))]).mean()
          
                          
# regions
print ('correlations regions')

regions = ['MidAtl','SouAtl','PacCon','PacNon','ENC','WNC','ESC','WSC','Mou','NewEng']

MidAtl = ['NY','NJ','PA','DE','MD','WA','VA','WV']
SouAtl = ['DE','MD','VA','WV','NC','SC','GA','FL','DC']
PacCon = ['CA','OR','WA']
PacNon = ['AK','HI']
ENC = ['IL','IN','MI','OH','WI']
WNC = ['IA','MN','MS','NE','ND','SD']
ESC = ['AL','KY','MS','TN']
WSC = ['AR','LA','OK','TX']
Mou = ['CO','WY','UT','NM','NV','ID','AZ','MT']
NewEng = ['CT','NH','ME','MA','RI','VT']

states_reg = [NewEng,MidAtl,ENC,WNC,SouAtl,ESC,WSC,Mou,PacCon,PacNon]

mc_reg_mer = []
mc_reg_era = []
for ind in range(len(states_reg)):
    wt_mer_reg = wt_mer[[state in states_reg[ind] for state in wt_mer.state]]
    wt_era_reg = wt_era[[state in states_reg[ind] for state in wt_era.state]]
    # unique values
    if len(wt_mer_reg.cor_id.unique())==1:
        mc_r_mer = 1
    else:
        mc_r_mer = np.array([list(c_mer[wt_mer_reg.cor_id.unique()[i],
                                        np.delete(wt_mer_reg.cor_id.unique(),i,0)]) for i in range(len(wt_mer_reg.cor_id.unique()))]).mean()
    if  len(wt_era_reg.cor_id.unique())==1:
        mc_r_era = 1
    else:
        mc_r_era = np.array([list(c_era[wt_era_reg.cor_id.unique()[i],
                                        np.delete(wt_era_reg.cor_id.unique(),i,0)]) for i in range(len(wt_era_reg.cor_id.unique()))]).mean()
    mc_reg_mer = mc_reg_mer + [mc_r_mer]
    mc_reg_era = mc_reg_era + [mc_r_era]


# states
print('correlations states')

mc_st_mer = []
mc_st_era = []

states = wt_mer.state.unique()
for state in states:
    wt_mer_s = wt_mer[wt_mer.state==state]
    wt_era_s = wt_era[wt_era.state==state]

    # unique values
    if len(wt_mer_s.cor_id.unique())==1:
        mc_s_mer = 1
    else:
        mc_s_mer = np.array([list(c_mer[wt_mer_s.cor_id.unique()[i],
                                        np.delete(wt_mer_s.cor_id.unique(),i,0)]) for i in range(len(wt_mer_s.cor_id.unique()))]).mean()
    if len(wt_era_s.cor_id.unique())==1:
        mc_s_era = 1
    else:
        mc_s_era = np.array([list(c_era[wt_era_s.cor_id.unique()[i],
                                        np.delete(wt_era_s.cor_id.unique(),i,0)]) for i in range(len(wt_era_s.cor_id.unique()))]).mean()
    mc_st_mer = mc_st_mer + [mc_s_mer]
    mc_st_era = mc_st_era + [mc_s_era]


# merge correlations
print('merge all correlations')
mean_cors_19y_wind = pd.DataFrame({'region': (['BPA','USA','NE'] + regions + states.tolist()) * 2,
                                   'dataset': ['MERRA2'] * 57 + ['ERA5'] * 57,
                                   'cor':[mc_BPA_mer,mc_USA_mer,mc_NE_mer] + mc_reg_mer + mc_st_mer +
                                          [mc_BPA_era,mc_USA_era,mc_NE_era] + mc_reg_era + mc_st_era})


# save wind speed correlations
print('save wind speed correlations')
mean_cors_19y_wind.to_csv(results_path + '/correlations_wind_19y.csv')