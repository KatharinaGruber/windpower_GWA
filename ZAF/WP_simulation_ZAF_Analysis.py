import numpy as np
import pandas as pd
import xarray as xr
from functools import reduce
from operator import mul
import sys
sys.path.append('./..')

from utils import stats
from paths_zaf import *


## Functions

def read_ZAFprod():
    '''
    function for reading production data from csv
     replace month names and convert to datetime format
    '''
    files = ['/Hourly_Distribution_data_NZ.csv',
             '/Hourly_Electricity_Production_[Load_Factor_[[%]]_data_EasternCape.csv',
             '/Hourly_Electricity_Production_[Load_Factor_[[%]]_data_NorthernCape.csv',
             '/Hourly_Electricity_Production_[Load_Factor_[[%]]_data_WesternCape.csv']
    regions = ['ZAF','Eastern Cape','Northern Cape','Western Cape']
    months = {' Januar ':'01.',
              ' Februar ':'02.',
              ' März ':'03.',
              ' April ':'04.',
              ' Mai ':'05.',
              ' Juni ':'06.',
              ' Juli ':'07.',
              ' August ':'08.',
              ' September ':'09.',
              ' Oktober ':'10.',
              ' November ':'11.',
              ' Dezember ':'12.'}
    prod_ZAF = []
    for (file,region) in zip(files,regions):
        ZAF = pd.read_csv(zaf_path + file,sep=';',parse_dates=[1])
        ZAF.columns =  ['year','date','technology','CF']
        ZAF = ZAF[ZAF.technology=='Onshore Wind\r\n'].drop(['technology','year'],axis=1)
        for month in months:
            ZAF.date = ZAF.date.str.replace(month,months.get(month)) # replace months by numbers
        prod_ZAF.append(pd.Series(ZAF.CF.values,
                                  index=pd.to_datetime(ZAF.date).values,
                                  name=region).sort_index())
    return(pd.concat(prod_ZAF,axis=1).tz_localize('Africa/Johannesburg'))

def load_results(dataset,gwa,regions):
    '''
    function for loading simulation results
     dataset is either MERRA2 or ERA5
     gwa is none, GWA2 or GWA3
    '''
    if gwa == 'GWA2':
        rpath = results_path + '/results_GWA2/'
    else:
        rpath = results_path + '/'
    if gwa == 'none':
        file = 'windpower_ZAF_'+dataset+'.nc'
    else:
        file = 'windpower_ZAF_'+dataset+'_GWA.nc'
    ZAF = xr.open_dataset(rpath + file).wp.to_dataframe().reset_index().set_index(['time','location']).unstack()
    # adapt datetime index of MERRA data (some hours are shifted by 30 min)
    if dataset == 'MERRA2':
        ZAF.index.values[ZAF.index.minute!=0] = ZAF.index.values[ZAF.index.minute!=0] - np.timedelta64(30,'m')
    # sum up per region
    ZAF = ZAF.groupby(regions,axis=1).sum(axis=1).tz_localize('UTC')
    ZAF['ZAF'] = ZAF.sum(axis=1)
    return(ZAF)

def get_cap_df(cap,comdate):
    '''
    function for getting hourly capacities
     cap is numpy array of installed capacities
     comdate is numpy array of commissioning dates in datetime format
    '''
    com = pd.DataFrame({'capacity': cap}).groupby(comdate).sum()
    cap_cum = com.capacity.cumsum()
    # if only years given for commissioning dates -> gradual capacity increase over year, full capacity at end of year
    dr = pd.date_range('1/1/2013','31/12/2019 23:00:00',freq = 'h')
    cap_ts = pd.Series(dr.map(cap_cum),index = dr)
    cap_ts[0] = cap_cum[cap_cum.index<=pd.Timestamp('2013-01-01')].max()
    if type(comdate[0]) == np.int64:
        return(cap_ts.interpolate(method='linear'))
    else:
        return(cap_ts.fillna(method='ffill'))

def gcdH(region):
    '''
    function to get capacity time series for a region
    '''
    cap = windparks[windparks.Area==region].Capacity.values
    com = windparks[windparks.Area==region].commissioning.astype(np.datetime64).values
    cdf = get_cap_df(cap,com).tz_localize('UTC')
    cdf.name = region
    return(cdf)

def analyse_ZAFh(region):
    '''
    analyse hourly wind power generation for a region
    '''
    comph = pd.DataFrame({'MERRA2':ZAFm[region],
                         'ERA5':ZAFe[region],
                         'MERRA2_GWA2':ZAFmg2[region],
                         'ERA5_GWA2':ZAFeg2[region],
                         'MERRA2_GWA3':ZAFmg3[region],
                         'ERA5_GWA3':ZAFeg3[region]})
    # get capacities
    caph = capdfH[region]
    # calculate capacity factors
    cfh = comph.div(caph,axis=0).tz_convert('Africa/Johannesburg')
    # add observed data
    cfh['obs'] = cfh.index.map(ZAFh[region])
    # remove capacity factors > 1 and lines with missing data
    cfh = cfh.mask(cfh>1).dropna()
    stat_h = pd.DataFrame({'ERA5':stats(cfh.ERA5,cfh.obs,False),
                           'ERA5_GWA2':stats(cfh.ERA5_GWA2,cfh.obs,False),
                           'ERA5_GWA3':stats(cfh.ERA5_GWA3,cfh.obs,False),
                           'MERRA2':stats(cfh.MERRA2,cfh.obs,False),
                           'MERRA2_GWA2':stats(cfh.MERRA2_GWA2,cfh.obs,False),
                           'MERRA2_GWA3':stats(cfh.MERRA2_GWA3,cfh.obs,False),
                           'obs':[np.nan,np.nan,np.nan,cfh.obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_h.columns = ['param','dataset',region]
    return(stat_h.set_index(['param','dataset']).transpose())

def analyse_ZAFd(region):
    '''
    analyse daily wind power generation for a region
    '''
    mask = (ZAFh[region].notna()*capdfH[region].notna()).replace(0,np.nan)
    comph = pd.DataFrame({'MERRA2':ZAFm[region].tz_convert('Africa/Johannesburg')*mask,
                         'ERA5':ZAFe[region].tz_convert('Africa/Johannesburg')*mask,
                         'MERRA2_GWA2':ZAFmg2[region].tz_convert('Africa/Johannesburg')*mask,
                         'ERA5_GWA2':ZAFeg2[region].tz_convert('Africa/Johannesburg')*mask,
                         'MERRA2_GWA3':ZAFmg3[region].tz_convert('Africa/Johannesburg')*mask,
                         'ERA5_GWA3':ZAFeg3[region].tz_convert('Africa/Johannesburg')*mask})
    # get capacities and mask
    caph = capdfH[region].tz_convert('Africa/Johannesburg')*mask
    # aggregate daily
    capd = caph.resample('D').sum()
    compd = comph.resample('D').sum()
    # calculate capacity factors
    cfd = compd.div(capd,axis=0)
    # add observed CFs
    cfd['obs'] = cfd.index.map((ZAFh[region]*mask).resample('D').mean())
    # remove capacity factors > 1 and missing data
    cfd = cfd.mask(cfd>1).dropna()
    stat_d = pd.DataFrame({'ERA5':stats(cfd.ERA5,cfd.obs,False),
                           'ERA5_GWA2':stats(cfd.ERA5_GWA2,cfd.obs,False),
                           'ERA5_GWA3':stats(cfd.ERA5_GWA3,cfd.obs,False),
                           'MERRA2':stats(cfd.MERRA2,cfd.obs,False),
                           'MERRA2_GWA2':stats(cfd.MERRA2_GWA2,cfd.obs,False),
                           'MERRA2_GWA3':stats(cfd.MERRA2_GWA3,cfd.obs,False),
                           'obs':[np.nan,np.nan,np.nan,cfd.obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_d.columns = ['param','dataset',region]
    return(stat_d.set_index(['param','dataset']).transpose())

def analyse_ZAFm(region):
    '''
    analyse monthly wind power generation for a region
    '''
    # mask for masking simulated data and capacities
    # (to only use timespans where also observed data are available)
    mask = (ZAFh[region].notna()*capdfH[region].notna()).replace(0,np.nan)
    comph = pd.DataFrame({'MERRA2':ZAFm[region].tz_convert('Africa/Johannesburg')*mask,
                         'ERA5':ZAFe[region].tz_convert('Africa/Johannesburg')*mask,
                         'MERRA2_GWA2':ZAFmg2[region].tz_convert('Africa/Johannesburg')*mask,
                         'ERA5_GWA2':ZAFeg2[region].tz_convert('Africa/Johannesburg')*mask,
                         'MERRA2_GWA3':ZAFmg3[region].tz_convert('Africa/Johannesburg')*mask,
                         'ERA5_GWA3':ZAFeg3[region].tz_convert('Africa/Johannesburg')*mask})
    # get capacities and mask
    caph = capdfH[region].tz_convert('Africa/Johannesburg')*mask
    # aggregate monthly
    capm = caph.resample('M').sum()
    compm = comph.resample('M').sum()
    # calculate capacity factors
    cfm = compm.div(capm,axis=0)
    # add observed data
    cfm['obs'] = cfm.index.map((ZAFh[region]*mask).resample('M').mean())
    # remove capacity factors > 1 and missing data
    cfm = cfm.mask(cfm>1).dropna()
    stat_m = pd.DataFrame({'ERA5':stats(cfm.ERA5,cfm.obs,False),
                           'ERA5_GWA2':stats(cfm.ERA5_GWA2,cfm.obs,False),
                           'ERA5_GWA3':stats(cfm.ERA5_GWA3,cfm.obs,False),
                           'MERRA2':stats(cfm.MERRA2,cfm.obs,False),
                           'MERRA2_GWA2':stats(cfm.MERRA2_GWA2,cfm.obs,False),
                           'MERRA2_GWA3':stats(cfm.MERRA2_GWA3,cfm.obs,False),
                           'obs':[np.nan,np.nan,np.nan,cfm.obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_m.columns = ['param','dataset',region]
    return(stat_m.set_index(['param','dataset']).transpose())

def tidy_res(results,temp):
    '''
    function for preparing results for merging
     results are results in wide format (one region per line)
     temp is temporal resolution (h,d,m)
    '''
    rt = results.transpose().reset_index().melt(id_vars=['param','dataset'], 
                                                value_vars=results.index, var_name='location')
    ts = pd.DataFrame({'temp':temp,'scale':'state'},
                      index = range(reduce(mul,results.shape)))
    ts.loc[rt.location=='ZAF','scale'] = 'country'
    return(pd.concat([rt,ts],axis=1))

## Analysis
# load windpark data
print('prepare windparks')
windparks = pd.read_csv(zaf_path + "/windparks_ZAF.csv", parse_dates=['commissioning'])

# load observed data
print('load observed data')
ZAFh = read_ZAFprod()

# get capacitiy time series for all parks
print('get capacities')
capes = ['Eastern Cape','Western Cape','Northern Cape']
capdfCH = pd.Series(capes,index = capes).apply(gcdH).transpose()
capdfZH = get_cap_df(windparks.Capacity.values,windparks.commissioning.astype(np.datetime64).values).tz_localize('UTC')
capdfZH.name = 'ZAF'
capdfH = pd.concat([capdfCH,capdfZH],axis = 1)

# load simulated data
print('load simulated data')
ZAFm = load_results('MERRA2','none',windparks.Area.values)
ZAFmg2 = load_results('MERRA2','GWA2',windparks.Area.values)
ZAFmg3 = load_results('MERRA2','GWA3',windparks.Area.values)
ZAFe = load_results('ERA5','none',windparks.Area.values)
ZAFeg2 = load_results('ERA5','GWA2',windparks.Area.values)
ZAFeg3 = load_results('ERA5','GWA3',windparks.Area.values)

# analyse results
# hourly
print('analyse hourly')
resZAFh = pd.concat(pd.Series(['ZAF']+capes).apply(analyse_ZAFh).to_list(),axis=0)
# daily
print('analyse daily')
resZAFd = pd.concat(pd.Series(np.unique(['ZAF']+capes)).apply(analyse_ZAFd).to_list(),axis=0)
# monthly
print('analyse monthly')
resZAFm = pd.concat(pd.Series(np.unique(['ZAF']+capes)).apply(analyse_ZAFm).to_list(),axis=0)

# tidy and merge results
print('tidy and merge results')
rZAFh = tidy_res(resZAFh,'h')
rZAFd = tidy_res(resZAFd,'d')
rZAFm = tidy_res(resZAFm,'m')
results = pd.concat([rZAFh,rZAFd,rZAFm],axis=0)
results['ds'] = results.dataset.str.split('_').apply(lambda x: x[0])
results['GWA'] = (results.dataset.str.split('_')+['none']).apply(lambda x: x[1])
# save results
print('save results')
results.to_csv(results_path + '/statZAF.csv')