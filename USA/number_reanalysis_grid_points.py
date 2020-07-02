#!/usr/bin/env python
# coding: utf-8

import glob
import numpy as np
import os
import pandas as pd
import xarray as xr
from scipy.stats import pearsonr

from paths_usa import *

from dask.diagnostics import ProgressBar
ProgressBar().register()

if results_path + '/number_grid_points.csv' not in glob.glob(results_path + '/*.csv'):
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


    # BPA
    print('BPA')

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

    # unique points BPA
    ngc_BPA_mer = len(BPA_mer.cor_id.unique())
    ngc_BPA_era = len(BPA_era.cor_id.unique())

    # USA
    print('USA')

    # unique points USA
    ngc_USA_mer = len(wt_mer.cor_id.unique())
    ngc_USA_era = len(wt_era.cor_id.unique())


    # New England
    print('New England')

    NE_states = ['CT','NH','ME','MA','RI','VT']
    wt_mer_NE = wt_mer[[state in NE_states for state in wt_mer.state]]
    wt_era_NE = wt_era[[state in NE_states for state in wt_era.state]]
    # unique mean correlations New England
    ngc_NE_mer = len(wt_mer_NE.cor_id.unique())
    ngc_NE_era = len(wt_era_NE.cor_id.unique())
              
                              
    # regions
    print ('regions')

    regions = ['MidAtl','SouAtl','PacCon','PacNon','ENC','WNC','ESC','WSC','Mou','NewEng']

    MidAtl = ['NY','NJ','PA']
    SouAtl = ['DE','MD','VA','WV','NC','SC','GA','FL','DC']
    PacCon = ['CA','OR','WA']
    PacNon = ['AK','HI']
    ENC = ['IL','IN','MI','OH','WI']
    WNC = ['IA','KS','MN','MS','NE','ND','SD']
    ESC = ['AL','KY','MS','TN']
    WSC = ['AR','LA','OK','TX']
    Mou = ['CO','WY','UT','NM','NV','ID','AZ','MT']
    NewEng = ['CT','NH','ME','MA','RI','VT']

    states_reg = [NewEng,MidAtl,ENC,WNC,SouAtl,ESC,WSC,Mou,PacCon,PacNon]

    ngc_reg_mer = []
    ngc_reg_era = []
    for ind in range(len(states_reg)):
        wt_mer_reg = wt_mer[[state in states_reg[ind] for state in wt_mer.state]]
        wt_era_reg = wt_era[[state in states_reg[ind] for state in wt_era.state]]
        # unique values
        ngc_r_mer = len(wt_mer_reg.cor_id.unique())
        ngc_r_era = len(wt_era_reg.cor_id.unique())
        ngc_reg_mer = ngc_reg_mer + [ngc_r_mer]
        ngc_reg_era = ngc_reg_era + [ngc_r_era]


    # states
    print('states')

    ngc_st_mer = []
    ngc_st_era = []

    states = wt_mer.state.unique()
    for state in states:
        wt_mer_s = wt_mer[wt_mer.state==state]
        wt_era_s = wt_era[wt_era.state==state]

        # unique values
        ngc_s_mer = len(wt_mer_s.cor_id.unique())
        ngc_s_era = len(wt_era_s.cor_id.unique())
        ngc_st_mer = ngc_st_mer + [ngc_s_mer]
        ngc_st_era = ngc_st_era + [ngc_s_era]


    # merge grid cell numbers
    print('merge all grid cell numbers')
    ngc_wind = pd.DataFrame({'region': (['BPA','USA','NE'] + regions + states.tolist()) * 2,
                                       'dataset': ['MERRA2'] * 57 + ['ERA5'] * 57,
                                       'cor':[ngc_BPA_mer,ngc_USA_mer,ngc_NE_mer] + ngc_reg_mer + ngc_st_mer +
                                              [ngc_BPA_era,ngc_USA_era,ngc_NE_era] + ngc_reg_era + ngc_st_era})


    # save
    print('save')
    ngc_wind.to_csv(results_path + '/number_grid_points.csv')