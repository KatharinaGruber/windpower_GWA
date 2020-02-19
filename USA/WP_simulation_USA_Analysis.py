#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import xarray as xr
import numpy as np
import datetime
import pickle as pkl
import pytz
import sys
sys.path.append('../')

from utils import stats

from dask.diagnostics import ProgressBar
ProgressBar().register()

from paths_usa import *


# define function for preparing capacity time series
def get_cap_df(cap,comdate):
    com = pd.DataFrame({'capacity': cap}).groupby(comdate).sum()
    cap_cum = com.capacity.cumsum()
    # if only years given for commissioning dates -> gradual capacity increase over year, full capacity at end of year
    if type(cap_cum.index.values[0]) == np.int64:
        cap_cum.index = [np.datetime64(str(int(year))+"-12-31 23:00:00") for year in cap_cum.index.values]
        # if missing years -> put capacity of year before
        drcc = pd.date_range(cap_cum.index[0],cap_cum.index[-1],freq = 'y')
        cap_cum = pd.Series(drcc.map(cap_cum),index = drcc).ffill()
    dr = pd.date_range('1/1/2000','31/12/2018 23:00:00',freq = 'h')
    cap_ts = pd.Series(dr.map(cap_cum),index = dr)
    cap_ts[0] = cap_cum[cap_cum.index<=pd.Timestamp('2000-01-01')].max()
    if type(comdate[0]) == np.int64:
        return(cap_ts.interpolate(method='linear'))
    else:
        return(cap_ts.fillna(method='ffill'))
    
    
### Prepare wind park data - Installed capacities

# get BPA windparks
BPA_parks = pd.read_csv(usa_path+"/BPA_windparks.csv")
# get windturbine locations/names
windturbines = pd.read_csv(usa_path+"/windturbines_usa.csv",delimiter=';')
wt = pd.read_csv(usa_path + "/turbine_data.csv",parse_dates = ['time'])
# get aggregated turbine data
turb_mer = pd.read_csv(usa_path + "/turbine_data_mer.csv",parse_dates = ['commissioning']).drop('Unnamed: 0',axis=1)
# find where BPA wind turbines are
pBPA = pd.DataFrame({'p': [park in BPA_parks.name.values for park in windturbines[windturbines.t_state!='GU'].p_name.values]})
# find where New England wind turbines are
NE_states = ['CT','NH','ME','MA','RI','VT']
NE_turbines = turb_mer[[state in ['CT','NH','ME','MA','RI','VT'] for state in turb_mer.state]]
# get capacities time series for all regions
cap_usa = get_cap_df(turb_mer.capacity.values,
                     pd.DatetimeIndex(turb_mer.commissioning).year.values).tz_localize('UTC').tz_convert('US/Central')
cap_IA = get_cap_df(turb_mer.capacity[turb_mer.state=='IA'].values,
                    pd.DatetimeIndex(turb_mer.commissioning[turb_mer.state=='IA']).year.values).tz_localize('UTC').tz_convert('US/Central')
cap_TX = get_cap_df(turb_mer.capacity[turb_mer.state=='TX'].values,
                    pd.DatetimeIndex(turb_mer.commissioning[turb_mer.state=='TX']).year.values).tz_localize('UTC').tz_convert('US/Central')
cap_BPA = get_cap_df(wt.capacity[pBPA.p].values,
                     pd.DatetimeIndex(wt.time[pBPA.p]).year.values).tz_localize('UTC')
cap_NE = get_cap_df(NE_turbines.capacity.values,
                    pd.DatetimeIndex(NE_turbines.commissioning).year.values).tz_localize('UTC').tz_convert('US/Eastern')
                    
# get IRENA capacities for observed generation USA
caps_irena = pd.read_csv(usa_path + '/IRENA_caps.csv').iloc[:,3:].T
caps_irena.columns = ['caps_MW']
# cumulative capacity
capc = caps_irena.caps_MW.str.replace(' ','').astype(np.int64).values
# added capacity
cap = np.append(capc[0],capc[1:]-capc[:-1])
comdate = caps_irena.index.values.astype(np.int64)
# get capacities time series
cap_usaIRENA = get_cap_df(cap,comdate).tz_localize('UTC').tz_convert('US/Central')

# aggregate daily or monthly where needed
cap_usam = cap_usa.resample('M').sum()
cap_IAm = cap_IA.resample('M').sum()
cap_TXm = cap_TX.resample('M').sum()
cap_TXd = cap_TX.resample('D').sum()
cap_TXh = cap_TX
cap_BPAm = cap_BPA.resample('M').sum()
cap_BPAd = cap_BPA.resample('D').sum()
cap_BPAh = cap_BPA
cap_NEm = cap_NE.resample('M').sum()
cap_usaIm = cap_usaIRENA.resample('M').sum()

### Analysis capacity factors

## USA monthly
# Load production data
# Source: https://www.eia.gov/electricity/data/browser/#/topic/0?agg=1,0,2&fuel=008&geo=vvvvvvvvvvvvo&sec=o3g&linechart=ELEC.GEN.WND-US-99.M~ELEC.GEN.WND-IA-99.M~ELEC.GEN.WND-TX-99.M&columnchart=ELEC.GEN.WND-US-99.M~ELEC.GEN.WND-IA-99.M~ELEC.GEN.WND-TX-99.M&map=ELEC.GEN.WND-US-99.M&freq=M&start=200101&end=201903&chartindexed=0&ctype=linechart&ltype=pin&rtype=s&pin=&rse=0&maptype=0
prod_USAm = pd.read_csv(usa_path+"/generation_data/USA_monthly/Net_generation_wind_all_sectors_monthly.csv",header=4)
# arrange data
# rename columns
prod_USAm.columns = ['time','wp_obs','Iowa','Texas']
# sort indices
prod_USAm = prod_USAm[~np.isnan(prod_USAm.wp_obs)].sort_index(axis=0 ,ascending=False)
# create datetime index
prod_USAm = prod_USAm.set_index(pd.to_datetime(prod_USAm.time.values)).drop(['time'],axis=1)
# cut after 2018
prod_USAm = prod_USAm[prod_USAm.index < np.datetime64("2019-01-01")].tz_localize('US/Central')
# Prepare simulated data
# load data
wpE = xr.open_dataset(results_path+"/windpower_USA_ERA5.nc").to_dataframe()
wpE_GWA = xr.open_dataset(results_path+"/windpower_USA_ERA5_GWA.nc").to_dataframe()
wpM = xr.open_dataset(results_path+"/windpower_USA_MERRA2.nc").to_dataframe()
wpM_GWA = xr.open_dataset(results_path+"/windpower_USA_MERRA2_GWA.nc").to_dataframe()
# merge data
wp_USA = pd.concat([wpE,wpE_GWA,wpM,wpM_GWA],axis=1).tz_localize('UTC').tz_convert('US/Central')
wp_USA.columns = ['ERA5','ERA5_GWA','MERRA2','MERRA2_GWA']
# aggregate monthly
wp_USAm = wp_USA.resample('M').sum()
# combine data and calculate capacity factors
cf_USAm = pd.concat([wp_USAm.div(cap_usam,axis=0),
                      prod_USAm.resample('M').sum().wp_obs*10**6/(cap_usaIm*10**3)],axis=1).dropna()[1:]
cf_USAm.columns = np.append(wp_USAm.columns,'wp_obs')
# Analyse
stats_USAm = pd.DataFrame({'ERA5':stats(cf_USAm.ERA5,cf_USAm.wp_obs,False),
                           'ERA5_GWA':stats(cf_USAm.ERA5_GWA,cf_USAm.wp_obs,False),
                           'MERRA2':stats(cf_USAm.MERRA2,cf_USAm.wp_obs,False),
                           'MERRA2_GWA':stats(cf_USAm.MERRA2_GWA,cf_USAm.wp_obs,False),
                           'obs':[np.nan,np.nan,np.nan,cf_USAm.wp_obs.mean()]},
                          index = ['cor','rmse','mbe','avg'])
stats_USAm_r = pd.DataFrame({'ERA5':stats(cf_USAm.ERA5,cf_USAm.wp_obs),
                             'ERA5_GWA':stats(cf_USAm.ERA5_GWA,cf_USAm.wp_obs),
                             'MERRA2':stats(cf_USAm.MERRA2,cf_USAm.wp_obs),
                             'MERRA2_GWA':stats(cf_USAm.MERRA2_GWA,cf_USAm.wp_obs),
                             'obs':[np.nan,np.nan,np.nan,round(cf_USAm.wp_obs.mean(),2)]},
                            index = ['cor','rmse','mbe','avg'])
# save statistical results
stats_USAm.to_csv(results_path+'/stats_USAm.csv')
stats_USAm_r.to_csv(results_path+'/stats_USAm_r.csv',sep=';')


## Iowa monthly
# Source: same as USA
# extract production of Iowa
prod_IAm = pd.DataFrame({'wp_obs':prod_USAm.Iowa},
                        index = prod_USAm.index)
# Prepare simulated data
# load data
wp_IAE = xr.open_dataset(results_path+"/windpower_states_ERA5.nc").sel({'state':"IA"}).to_dataframe().drop('state',axis=1)
wp_IAE_GWA = xr.open_dataset(results_path+"/windpower_states_ERA5_GWA.nc").sel({'state':"IA"}).to_dataframe().drop('state',axis=1)
wp_IAM = xr.open_dataset(results_path+"/windpower_states_MERRA2.nc").sel({'state':"IA"}).to_dataframe().drop('state',axis=1)
wp_IAM_GWA = xr.open_dataset(results_path+"/windpower_states_MERRA2_GWA.nc").sel({'state':"IA"}).to_dataframe().drop('state',axis=1)
# merge data
wp_IAh = pd.concat([wp_IAE,wp_IAE_GWA,wp_IAM,wp_IAM_GWA],axis=1).tz_localize('UTC').tz_convert('US/Central')
wp_IAh.columns = ['ERA5','ERA5_GWA','MERRA2','MERRA2_GWA']
# sum up monthly
comp_IAm = wp_IAh.resample('M').sum()
# add production data
comp_IAm = pd.concat([comp_IAm,prod_IAm.resample('M').sum().wp_obs*10**6],axis=1)
# calculate capacity factors
cf_IAm = comp_IAm.div(cap_IAm,axis=0).dropna()
# Analyse
stats_IAm = pd.DataFrame({'ERA5':stats(cf_IAm.ERA5,cf_IAm.wp_obs,False),
                          'ERA5_GWA':stats(cf_IAm.ERA5_GWA,cf_IAm.wp_obs,False),
                          'MERRA2':stats(cf_IAm.MERRA2,cf_IAm.wp_obs,False),
                          'MERRA2_GWA':stats(cf_IAm.MERRA2_GWA,cf_IAm.wp_obs,False),
                          'obs':[np.nan,np.nan,np.nan,cf_IAm.wp_obs.mean()]},
                         index = ['cor','rmse','mbe','avg'])
stats_IAm_r = pd.DataFrame({'ERA5':stats(cf_IAm.ERA5,cf_IAm.wp_obs),
                            'ERA5_GWA':stats(cf_IAm.ERA5_GWA,cf_IAm.wp_obs),
                            'MERRA2':stats(cf_IAm.MERRA2,cf_IAm.wp_obs),
                            'MERRA2_GWA':stats(cf_IAm.MERRA2_GWA,cf_IAm.wp_obs),
                            'obs':[np.nan,np.nan,np.nan,round(cf_IAm.wp_obs.mean(),2)]},
                           index = ['cor','rmse','mbe','avg'])
# save statistical results
stats_IAm.to_csv(results_path+'/stats_IAm.csv')
stats_IAm_r.to_csv(results_path+'/stats_IAm_r.csv',sep=';')


## Texas monthly
# Source: same as USA
# extract production of Texas
prod_TXm = pd.DataFrame({'wp_obs':prod_USAm.Texas},
                        index = prod_USAm.index)
# Prepare simulated data
# load data
wp_TXE = xr.open_dataset(results_path+"/windpower_states_ERA5.nc").sel({'state':"TX"}).to_dataframe().drop('state',axis=1)
wp_TXE_GWA = xr.open_dataset(results_path+"/windpower_states_ERA5_GWA.nc").sel({'state':"TX"}).to_dataframe().drop('state',axis=1)
wp_TXM = xr.open_dataset(results_path+"/windpower_states_MERRA2.nc").sel({'state':"TX"}).to_dataframe().drop('state',axis=1)
wp_TXM_GWA = xr.open_dataset(results_path+"/windpower_states_MERRA2_GWA.nc").sel({'state':"TX"}).to_dataframe().drop('state',axis=1)
# merge data
wp_TXh = pd.concat([wp_TXE,wp_TXE_GWA,wp_TXM,wp_TXM_GWA],axis=1).tz_localize('UTC').tz_convert('US/Central')
wp_TXh.columns = ['ERA5','ERA5_GWA','MERRA2','MERRA2_GWA']
# aggregate monthly
comp_TXm = wp_TXh.resample('M').sum()
# add production data
comp_TXm = pd.concat([comp_TXm,prod_TXm.resample('M').sum().wp_obs*10**6],axis=1)
# calculate capacity factors
cf_TXm = comp_TXm.div(cap_TXm,axis=0).dropna()
# Analyse
stats_TXm = pd.DataFrame({'ERA5':stats(cf_TXm.ERA5,cf_TXm.wp_obs,False),
                          'ERA5_GWA':stats(cf_TXm.ERA5_GWA,cf_TXm.wp_obs,False),
                          'MERRA2':stats(cf_TXm.MERRA2,cf_TXm.wp_obs,False),
                          'MERRA2_GWA':stats(cf_TXm.MERRA2_GWA,cf_TXm.wp_obs,False),
                          'obs':[np.nan,np.nan,np.nan,cf_TXm.wp_obs.mean()]},
                         index = ['cor','rmse','mbe','avg'])
stats_TXm_r = pd.DataFrame({'ERA5':stats(cf_TXm.ERA5,cf_TXm.wp_obs),
                          'ERA5_GWA':stats(cf_TXm.ERA5_GWA,cf_TXm.wp_obs),
                          'MERRA2':stats(cf_TXm.MERRA2,cf_TXm.wp_obs),
                          'MERRA2_GWA':stats(cf_TXm.MERRA2_GWA,cf_TXm.wp_obs),
                          'obs':[np.nan,np.nan,np.nan,round(cf_TXm.wp_obs.mean(),2)]},
                         index = ['cor','rmse','mbe','avg'])
# save statistical results
stats_TXm.to_csv(results_path+'/stats_TXm.csv')
stats_TXm_r.to_csv(results_path+'/stats_TXm_r.csv',sep=';')


## Texas hourly
# Load production data
# Source: http://mis.ercot.com/misapp/GetReports.do?reportTypeId=13424&reportTitle=Hourly%20Aggregated%20Wind%20Output&showHTMLView=&mimicKey
files =  os.listdir(usa_path+"/generation_data/TX_hourly")
for file in files:
    if 'prod_TXh_xl' in globals():
        prod_TXh_xl = pd.concat([prod_TXh_xl,pd.read_excel(usa_path+"/generation_data/TX_hourly/"+file,sheet_name="numbers")],axis=0)
    else:
        prod_TXh_xl = pd.read_excel(usa_path+"/generation_data/TX_hourly/"+file,sheet_name="numbers")
# arrange data - get production data
# and consider daylight saving time and convert timezone to UTC
# Prepare simulated data
prod_TXh = pd.DataFrame({'wp_obs':prod_TXh_xl['Total Wind Output, MW'].values},
                         index = pd.to_datetime(prod_TXh_xl.Date.values)).tz_localize('US/Central',ambiguous='infer').tz_convert('UTC')
# remove duplicates
prod_TXh = prod_TXh[~prod_TXh.index.duplicated()]*10**3
# load simulated data
wp_TXE = xr.open_dataset(results_path+"/windpower_states_ERA5.nc").sel(state='TX').drop('state').to_dataframe()
wp_TXE_GWA = xr.open_dataset(results_path+"/windpower_states_ERA5_GWA.nc").sel(state='TX').drop('state').to_dataframe()
wp_TXM = xr.open_dataset(results_path+"/windpower_states_MERRA2.nc").sel(state='TX').drop('state').to_dataframe()
wp_TXM_GWA = xr.open_dataset(results_path+"/windpower_states_MERRA2_GWA.nc").sel(state='TX').drop('state').to_dataframe()
# merge
comp_TXh = pd.concat([wp_TXE,wp_TXE_GWA,wp_TXM,wp_TXM_GWA],axis=1)
comp_TXh.columns = ['ERA5','ERA5_GWA','MERRA2','MERRA2_GWA']
comp_TXh = comp_TXh.tz_localize('UTC')
# add production data
comp_TXh['wp_obs'] = comp_TXh.index.map(prod_TXh.wp_obs)
comp_TXh = comp_TXh.fillna(method='ffill').dropna() # remove rows before 2016 and one missing value (2016-07-07 17:00:00)
# calculate capacity factor
cf_TXh = comp_TXh.div(cap_TXh.tz_convert('UTC'),axis=0).dropna()
# Analyse
stats_TXh = pd.DataFrame({'ERA5':stats(cf_TXh.ERA5,cf_TXh.wp_obs,False),
                          'ERA5_GWA':stats(cf_TXh.ERA5_GWA,cf_TXh.wp_obs,False),
                          'MERRA2':stats(cf_TXh.MERRA2,cf_TXh.wp_obs,False),
                          'MERRA2_GWA':stats(cf_TXh.MERRA2_GWA,cf_TXh.wp_obs,False),
                          'obs':[np.nan,np.nan,np.nan,cf_TXh.wp_obs.mean()]},
                         index = ['cor','rmse','mbe','avg'])
stats_TXh_r = pd.DataFrame({'ERA5':stats(cf_TXh.ERA5,cf_TXh.wp_obs),
                            'ERA5_GWA':stats(cf_TXh.ERA5_GWA,cf_TXh.wp_obs),
                            'MERRA2':stats(cf_TXh.MERRA2,cf_TXh.wp_obs),
                            'MERRA2_GWA':stats(cf_TXh.MERRA2_GWA,cf_TXh.wp_obs),
                            'obs':[np.nan,np.nan,np.nan,round(cf_TXh.wp_obs.mean(),2)]},
                           index = ['cor','rmse','mbe','avg'])
# save statistical results
stats_TXh.to_csv(results_path+'/stats_TXh.csv')
stats_TXh_r.to_csv(results_path+'/stats_TXh_r.csv',sep=';')


## Texas Daily
# Prepare data
# aggregate per day
comp_TXd = comp_TXh.tz_convert('US/Central').resample('D').sum()
# calculate capacity factors
cf_TXd = comp_TXd.div(cap_TXd,axis=0).dropna()
# Analyse
stats_TXd = pd.DataFrame({'ERA5':stats(cf_TXd.ERA5,cf_TXd.wp_obs,False),
                          'ERA5_GWA':stats(cf_TXd.ERA5_GWA,cf_TXd.wp_obs,False),
                          'MERRA2':stats(cf_TXd.MERRA2,cf_TXd.wp_obs,False),
                          'MERRA2_GWA':stats(cf_TXd.MERRA2_GWA,cf_TXd.wp_obs,False),
                          'obs':[np.nan,np.nan,np.nan,cf_TXd.wp_obs.mean()]},
                         index = ['cor','rmse','mbe','avg'])
stats_TXd_r = pd.DataFrame({'ERA5':stats(cf_TXd.ERA5,cf_TXd.wp_obs),
                            'ERA5_GWA':stats(cf_TXd.ERA5_GWA,cf_TXd.wp_obs),
                            'MERRA2':stats(cf_TXd.MERRA2,cf_TXd.wp_obs),
                            'MERRA2_GWA':stats(cf_TXd.MERRA2_GWA,cf_TXd.wp_obs),
                            'obs':[np.nan,np.nan,np.nan,round(cf_TXd.wp_obs.mean(),2)]},
                           index = ['cor','rmse','mbe','avg'])
# save statistical results
stats_TXd.to_csv(results_path+'/stats_TXd.csv')
stats_TXd_r.to_csv(results_path+'/stats_TXd_r.csv',sep=';')


## BPA hourly + daily + monthly
# Load production data
# Source: https://transmission.bpa.gov/Business/Operations/Wind/  
# arrange data - get production data
# and consider daylight saving time and convert timezone to UTC
# load and prepare production data in convenient format
# for years 2007 to 2010 files are not formatted equally, so define columns to select
col = pd.DataFrame({'dat1':[2,2,3,3],
                    'dat2':[10,11,12,12],
                    'time1':[1,1,1,1],
                    'time2':[9,10,10,10],
                    'lines':[14,16,16,16]},
                   index = ['07','08','09','10'])
# for years 2011 to 2019 only lines are different
lines = pd.DataFrame({'x': [18,18,18,18,18,18,21,19,19]},
                     index = [11,12,13,14,15,16,17,18,19])
for year in ['07','08','09']:
    BPAxl = pd.read_excel(usa_path+"/generation_data/BPA_5min/TotalWindLoad_5Min_"+year+".xls",
                          skiprows=col.lines[year])
    # for first half of year less rows - drop them
    BPAxl_drop = BPAxl.dropna(subset=[BPAxl.columns[1]],axis=0)
    # merge data for first and second half of year in datetime indexed data frame
    BPAxl_dat = pd.DataFrame({'wp_obs':pd.concat([BPAxl_drop.iloc[:,col.dat1[year]],
                                                  BPAxl.iloc[:,col.dat2[year]]]).values},
                             index = pd.concat([pd.to_datetime(BPAxl_drop.iloc[:,col.time1[year]]),
                                                pd.to_datetime(BPAxl.iloc[:,col.time2[year]])]).values)
    # there seems to be some numeric error, therefore round index to minutes
    BPAxl_dat.index = BPAxl_dat.index.round("Min")
    
    # add the timezone and interpolate missing data (a few minutes are missing)
    #   create a dataframe with 5 minute timesteps to match values
    #   this is necessary, as there is lines for hours omitted due to daylight saving
    #   but no additional hours due to duplication of an hour in daylight saving in the second half of the year
    #   therefore the ambiguous argument 'infer' does not work
    #   by mapping the values of this hour are duplicated and then infer works
    ind = pd.date_range(start='1/1/20'+year+' 00:00:00',
                        end='31/12/20'+year+' 23:59:00', freq='5min',tz="US/Pacific")
    ind2 = [i.replace(tzinfo=None) for i in ind]
    BPAdat = pd.DataFrame({'wp_obs':np.nan},
                          index = ind2)
    BPAdat.wp_obs = BPAdat.index.map(BPAxl_dat.wp_obs)
    BPAtz = BPAdat.tz_localize('US/Pacific',ambiguous='infer').interpolate()
    # concatenate years
    if 'prod_BPA5min' in globals():
        prod_BPA5min = pd.concat([prod_BPA5min,BPAtz],axis=0)
    else:
        prod_BPA5min = BPAtz
# year 2010 is the only of the first set of years that actually has data for the duplicate hour
# therefore needs to be handled differently
year = '10'
BPAxl = pd.read_excel(usa_path+"/generation_data/BPA_5min/TotalWindLoad_5Min_"+year+".xls",
                      skiprows=col.lines[year])
BPAxl_drop = BPAxl.dropna(subset=[BPAxl.columns[1]],axis=0)
BPAxl_dat = pd.DataFrame({'wp_obs':pd.concat([BPAxl_drop.iloc[:,col.dat1[year]],
                                              BPAxl.iloc[:,col.dat2[year]]]).values},
                         index = pd.concat([pd.to_datetime(BPAxl_drop.iloc[:,col.time1[year]]),
                                            pd.to_datetime(BPAxl.iloc[:,col.time2[year]])]).values)
BPAtz = BPAxl_dat.tz_localize('US/Pacific',ambiguous='infer').interpolate()
prod_BPA5min = pd.concat([prod_BPA5min,BPAtz],axis=0)
for year in range(11,19):
    # read data from two sheets separately
    BPAxl1 = pd.read_excel(usa_path+"/generation_data/BPA_5min/WindGenTotalLoadYTD_20"+str(year)+".xls",
                           sheet_name="January-June",
                           skiprows=lines.x[year])
    BPAxl2 = pd.read_excel(usa_path+"/generation_data/BPA_5min/WindGenTotalLoadYTD_20"+str(year)+".xls",
                           sheet_name="July-December",
                           skiprows=lines.x[year])
    # merge data for first and second half of year in datetime indexed data frame
    # and add a timezone and interpolate missing values
    BPAdat = pd.DataFrame({'wp_obs':pd.concat([BPAxl1.iloc[:,2],BPAxl2.iloc[:,2]]).values},
                          index = pd.concat([pd.to_datetime(BPAxl1.iloc[:,0]),
                                             pd.to_datetime(BPAxl2.iloc[:,0])]).values
                         ).tz_localize('US/Pacific',ambiguous='infer').interpolate()
    # concatenate years
    prod_BPA5min = pd.concat([prod_BPA5min,BPAdat],axis=0)
# select only data of interest and convert to UTC
prod_BPA5min = prod_BPA5min[prod_BPA5min.index<pd.to_datetime("2019-01-01",utc=True)].tz_convert('UTC')
# aggregate hourly  
prod_BPAh = prod_BPA5min.resample('H').sum()/12* 10**3

# Prepare simulated data
# load data
wp_BPAE = xr.open_dataset(results_path+"/windpower_BPA_ERA5.nc").to_dataframe()
wp_BPAE_GWA = xr.open_dataset(results_path+"/windpower_BPA_ERA5_GWA.nc").to_dataframe()
wp_BPAM = xr.open_dataset(results_path+"/windpower_BPA_MERRA2.nc").to_dataframe()
wp_BPAM_GWA = xr.open_dataset(results_path+"/windpower_BPA_MERRA2_GWA.nc").to_dataframe()
# merge data
comp_BPAh = pd.concat([wp_BPAE,wp_BPAE_GWA,wp_BPAM,wp_BPAM_GWA],axis=1)
comp_BPAh.columns = ['ERA5','ERA5_GWA','MERRA2','MERRA2_GWA']
# set timezone
comp_BPAh = comp_BPAh.tz_localize('UTC')
# add production data
comp_BPAh['wp_obs'] = comp_BPAh.index.map(prod_BPAh.wp_obs)
# remove lines before production
comp_BPAh = comp_BPAh.dropna()


## BPA hourly
cf_BPAh = comp_BPAh.div(cap_BPAh,axis=0).dropna()
# Analyse
stats_BPAh = pd.DataFrame({'ERA5':stats(cf_BPAh.ERA5,cf_BPAh.wp_obs,False),
                          'ERA5_GWA':stats(cf_BPAh.ERA5_GWA,cf_BPAh.wp_obs,False),
                          'MERRA2':stats(cf_BPAh.MERRA2,cf_BPAh.wp_obs,False),
                          'MERRA2_GWA':stats(cf_BPAh.MERRA2_GWA,cf_BPAh.wp_obs,False),
                          'obs':[np.nan,np.nan,np.nan,cf_BPAh.wp_obs.mean()]},
                         index = ['cor','rmse','mbe','avg'])
stats_BPAh_r = pd.DataFrame({'ERA5':stats(cf_BPAh.ERA5,cf_BPAh.wp_obs),
                            'ERA5_GWA':stats(cf_BPAh.ERA5_GWA,cf_BPAh.wp_obs),
                            'MERRA2':stats(cf_BPAh.MERRA2,cf_BPAh.wp_obs),
                            'MERRA2_GWA':stats(cf_BPAh.MERRA2_GWA,cf_BPAh.wp_obs),
                            'obs':[np.nan,np.nan,np.nan,round(cf_BPAh.wp_obs.mean(),2)]},
                           index = ['cor','rmse','mbe','avg'])
# save statistical results
stats_BPAh.to_csv(results_path+'/stats_BPAh.csv')
stats_BPAh_r.to_csv(results_path+'/stats_BPAh_r.csv',sep=';')

## BPA daily
# aggregate daily
comp_BPAd = comp_BPAh.resample('D').sum()[1:] # drop first day because incomplete
# calculate capacity factors
cf_BPAd = comp_BPAd.div(cap_BPAd,axis=0).dropna()
# Analyse
stats_BPAd = pd.DataFrame({'ERA5':stats(cf_BPAd.ERA5,cf_BPAd.wp_obs,False),
                          'ERA5_GWA':stats(cf_BPAd.ERA5_GWA,cf_BPAd.wp_obs,False),
                          'MERRA2':stats(cf_BPAd.MERRA2,cf_BPAd.wp_obs,False),
                          'MERRA2_GWA':stats(cf_BPAd.MERRA2_GWA,cf_BPAd.wp_obs,False),
                          'obs':[np.nan,np.nan,np.nan,cf_BPAd.wp_obs.mean()]},
                         index = ['cor','rmse','mbe','avg'])
stats_BPAd_r = pd.DataFrame({'ERA5':stats(cf_BPAd.ERA5,cf_BPAd.wp_obs),
                            'ERA5_GWA':stats(cf_BPAd.ERA5_GWA,cf_BPAd.wp_obs),
                            'MERRA2':stats(cf_BPAd.MERRA2,cf_BPAd.wp_obs),
                            'MERRA2_GWA':stats(cf_BPAd.MERRA2_GWA,cf_BPAd.wp_obs),
                            'obs':[np.nan,np.nan,np.nan,round(cf_BPAd.wp_obs.mean(),2)]},
                           index = ['cor','rmse','mbe','avg'])
# save statistical results
stats_BPAd.to_csv(results_path+'/stats_BPAd.csv')
stats_BPAd_r.to_csv(results_path+'/stats_BPAd_r.csv',sep=';')

## BPA monthly
# aggregate monthly
comp_BPAm = comp_BPAh.resample('M').sum()
# calculate capacity factors
cf_BPAm = comp_BPAm.div(cap_BPAm,axis=0).dropna()
# Analyse
stats_BPAm = pd.DataFrame({'ERA5':stats(cf_BPAm.ERA5,cf_BPAm.wp_obs,False),
                          'ERA5_GWA':stats(cf_BPAm.ERA5_GWA,cf_BPAm.wp_obs,False),
                          'MERRA2':stats(cf_BPAm.MERRA2,cf_BPAm.wp_obs,False),
                          'MERRA2_GWA':stats(cf_BPAm.MERRA2_GWA,cf_BPAm.wp_obs,False),
                          'obs':[np.nan,np.nan,np.nan,cf_BPAm.wp_obs.mean()]},
                         index = ['cor','rmse','mbe','avg'])
stats_BPAm_r = pd.DataFrame({'ERA5':stats(cf_BPAm.ERA5,cf_BPAm.wp_obs),
                            'ERA5_GWA':stats(cf_BPAm.ERA5_GWA,cf_BPAm.wp_obs),
                            'MERRA2':stats(cf_BPAm.MERRA2,cf_BPAm.wp_obs),
                            'MERRA2_GWA':stats(cf_BPAm.MERRA2_GWA,cf_BPAm.wp_obs),
                            'obs':[np.nan,np.nan,np.nan,round(cf_BPAm.wp_obs.mean(),2)]},
                           index = ['cor','rmse','mbe','avg'])
# save statistical results
stats_BPAm.to_csv(results_path+'/stats_BPAm.csv')
stats_BPAm_r.to_csv(results_path+'/stats_BPAm_r.csv',sep=';')


## New England monthly
# Load production data
# Source: https://www.iso-ne.com/isoexpress/web/reports/load-and-demand/-/tree/sys-load-eei-fmt
# years 2000 - 2015 in one file
for year in range(2000,2016):
    NExl = pd.read_excel(usa_path+"/generation_data/NewEngland_monthly/2000-2015_energy_peak_source.xls",
                         sheet_name = str(year))
    windy = NExl.iloc[np.where(NExl.iloc[:,0]=='Wind')[0][0]][2:14].values
    if 'windally' in globals():
        windally = windally + list(windy)
    else:
        windally = list(windy)
for year in range(2016,2019):
    NExl = pd.read_excel(usa_path+"/generation_data/NewEngland_monthly/"+str(year)+"_energy_peak_by_source.xlsx",
                         sheet_name = str(year))
    windy = NExl.iloc[np.where(NExl.iloc[:,0]=='WIND')[0][0]][1:13].values
    windally = windally + list(windy)
prod_NE = pd.DataFrame({'wp_obs':windally},
                       index = pd.date_range(start='01-01-2000',
                                             end='01-01-2019',
                                             freq='m').values.astype('datetime64[M]').astype('datetime64[s]'))# type of index is changed to get beginning of month
prod_NE = prod_NE.tz_localize('US/Eastern')
# Prepare simulated data
# load data
wp_NEE = xr.open_dataset(results_path+"/windpower_NewEngland_ERA5.nc").to_dataframe()
wp_NEE_GWA = xr.open_dataset(results_path+"/windpower_NewEngland_ERA5_GWA.nc").to_dataframe()
wp_NEM = xr.open_dataset(results_path+"/windpower_NewEngland_MERRA2.nc").to_dataframe()
wp_NEM_GWA = xr.open_dataset(results_path+"/windpower_NewEngland_MERRA2_GWA.nc").to_dataframe()
# merge data
comp_NEh = pd.concat([wp_NEE,wp_NEE_GWA,wp_NEM,wp_NEM_GWA],axis=1)
comp_NEh.columns = ['ERA5','ERA5_GWA','MERRA2','MERRA2_GWA']
comp_NEm = comp_NEh.tz_localize('UTC').tz_convert('US/Eastern').resample('M').sum()[1:-1] # cut last and first because incomplete
# add production data
comp_NEm['wp_obs'] = comp_NEm.index.map(prod_NE[prod_NE.wp_obs>0].resample('M').mean().wp_obs) * 10**6
comp_NEm = comp_NEm.dropna()
# calculate capacity factor
cf_NEm = comp_NEm.div(comp_NEm.index.map(cap_NEm),axis=0).dropna()
# Analyse
stats_NEm = pd.DataFrame({'ERA5':stats(cf_NEm.ERA5,cf_NEm.wp_obs,False),
                          'ERA5_GWA':stats(cf_NEm.ERA5_GWA,cf_NEm.wp_obs,False),
                          'MERRA2':stats(cf_NEm.MERRA2,cf_NEm.wp_obs,False),
                          'MERRA2_GWA':stats(cf_NEm.MERRA2_GWA,cf_NEm.wp_obs,False),
                          'obs':[np.nan,np.nan,np.nan,cf_NEm.wp_obs.mean()]},
                         index = ['cor','rmse','mbe','avg'])
stats_NEm_r = pd.DataFrame({'ERA5':stats(cf_NEm.ERA5,cf_NEm.wp_obs),
                            'ERA5_GWA':stats(cf_NEm.ERA5_GWA,cf_NEm.wp_obs),
                            'MERRA2':stats(cf_NEm.MERRA2,cf_NEm.wp_obs),
                            'MERRA2_GWA':stats(cf_NEm.MERRA2_GWA,cf_NEm.wp_obs),
                            'obs':[np.nan,np.nan,np.nan,round(cf_NEm.wp_obs.mean(),2)]},
                           index = ['cor','rmse','mbe','avg'])
# save statistical results
stats_NEm.to_csv(results_path+'/stats_NEm.csv')
stats_NEm_r.to_csv(results_path+'/stats_NEm_r.csv',sep=';')
