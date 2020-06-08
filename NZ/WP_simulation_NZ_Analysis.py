import numpy as np
import pandas as pd
import xarray as xr
from functools import reduce
from operator import mul
import sys
sys.path.append('./..')

from utils import stats
from paths_nz import *


## Functions

def read_NZprod():
    '''
    function for reading production data from csv
     read files per month separately and join
    '''
    years = range(1997,2020)
    months = range(1,13)
    for year in years:
        for month in months:
            try:
                prod_NZm = pd.read_csv(nz_path + "/generation/"+str(year)+f"{month:02d}"+"_Generation_MD.csv")
            except:
                continue
            if("prod_NZ" in locals()):
                prod_NZ = pd.concat([prod_NZ,prod_NZm[prod_NZm.Fuel_Code=="Wind"]])
            else:
                prod_NZ = prod_NZm[prod_NZm.Fuel_Code=="Wind"]
    return(prod_NZ)

def prep_gen(loc):
    '''
    helper function for preparing productiond data
     selects data of one wind park loc
     transforms the trading periods to datetime format
     and aggregates data hourly
    '''
    prod_loc = prod_NZts[prod_NZts.POC_Gen==loc].drop(['POC_Gen','POC_Code'],axis=1)
    ind = pd.date_range(start=pd.to_datetime(prod_loc.Trading_date.values[0]),
                        freq='H',
                        periods=len(prod_loc.TP)/2,
                        tz='NZ').tz_convert('UTC').repeat(2)
    prod_loc = prod_loc.set_index(ind)
    prod_loch = prod_loc.prod_kW.resample('H').sum()
    return(prod_loch)

def rm_constTS(wpt,lim=24):
    '''
    function for removing constant parts of time series
     all series of more than lim (standard: 24 (hours))
     are removed from the dataset
    '''
    wpt1 = wpt.copy(deep=True)
    wpt1.index = wpt.index - np.timedelta64(1,'h')
    # starts of constant timeseries
    s = np.where((((wpt-wpt1).values[1:]==0).astype(int)-
                 ((wpt-wpt1).values[:-1]==0).astype(int))==1)[0]
    # ends of constant timeseries
    e = np.where((((wpt-wpt1).values[1:]==0).astype(int)-
                 ((wpt-wpt1).values[:-1]==0).astype(int))==-1)[0]
    # filter starts and ends of rows of constant that are longer than 24 hours
    sd = s[np.where((e-s)>lim)[0]]
    ed = e[np.where((e-s)>lim)[0]]
    rmdf = pd.Series(0,index=wpt.index)
    for i in range(len(sd)):
        rmdf.iloc[sd[i]:ed[i]] = 1
    return(wpt.where(rmdf==0))

def tidy_prod(prod_NZ):
    # bring from wide to tidy format
    prod_NZt = pd.melt(prod_NZ,
                       id_vars = ['POC_Code','Gen_Code','Trading_date'],
                       value_vars = ['TP' + str(i) for i in range(1,51)],
                       value_name="prod_kW").dropna()

    # extract number of trading period and sort
    prod_NZt['TP'] = [int(tp[2:]) for tp in prod_NZt.variable.values]
    prod_NZts = prod_NZt.sort_values(by=['POC_Code','Gen_Code','Trading_date','TP'])           
    # separate data by datetimeindex and location
    # as neither POC Codes (te rere hau & twf 3) nor Gen Codes (twf 12) are unique - combine to get unique values per location
    prod_NZts['POC_Gen'] = prod_NZts.Gen_Code + prod_NZts.POC_Code 
    prod_NZ = pd.Series(prod_NZts.POC_Gen.unique()).apply(prep_gen).transpose()
    prod_NZ.columns = prod_NZts.POC_Gen.unique()
    # sum up both west wind
    prod_NZ['west_wind'] = prod_NZ.west_windWWD1102 + prod_NZ.west_windWWD1103
    prod_NZ = prod_NZ.drop(['west_windWWD1102','west_windWWD1103'],axis=1)
    # associate wind park names
    # prod_NZ[1] probaly is twf 2 because it has mostly higher production than prod_NZ[0], which has lower installed capacity
    gen_wps = [n[:-7] for n in prod_NZ.columns]
    gen_wps[0] = "twf_1"
    gen_wps[1] = "twf_2"
    gen_wps[7] = "west_wind"
    prod_NZ.columns = gen_wps
    # remove constant parts of time series longer than one day
    prod_NZh = prod_NZ.apply(rm_constTS,axis=0)
    return(prod_NZh)

def load_results(dataset,gwa,parks):
    '''
    function for loading simulation results
     dataset is either MERRA2 or ERA5
     gwa is none, GWA2 or GWA3
     parks are the parkname codes
    '''
    if gwa == 'GWA2':
        rpath = results_path + '/results_GWA2/'
    else:
        rpath = results_path + '/'
    if gwa == 'none':
        file = 'windpower_NZ_'+dataset+'.nc'
    else:
        file = 'windpower_NZ_'+dataset+'_GWA.nc'
    NZ = xr.open_dataset(rpath + file).wp.to_dataframe().reset_index().set_index(['time','location']).unstack()
    # adapt datetime index of MERRA data (some hours are shifted by 30 min)
    if dataset == 'MERRA2':
        NZ.index.values[NZ.index.minute!=0] = NZ.index.values[NZ.index.minute!=0] - np.timedelta64(30,'m')
    # sum up Te Rere Hau and adapt names to generation data
    NZ = NZ.groupby(parks,axis=1).sum(axis=1).tz_localize('UTC')
    return(NZ)

def get_cap_df(cap,comdate):
    '''
    function for getting hourly capacities
     cap is numpy array of installed capacities
     comdate is numpy array of commissioning dates in datetime format
    '''
    com = pd.DataFrame({'capacity': cap}).groupby(comdate).sum()
    cap_cum = com.capacity.cumsum()
    # if only years given for commissioning dates -> gradual capacity increase over year, full capacity at end of year
    dr = pd.date_range('1/1/1997','31/12/2019 23:00:00',freq = 'h')
    cap_ts = pd.Series(dr.map(cap_cum),index = dr)
    cap_ts[0] = cap_cum[cap_cum.index<=pd.Timestamp('1997-01-01')].max()
    if type(comdate[0]) == np.int64:
        return(cap_ts.interpolate(method='linear'))
    else:
        return(cap_ts.fillna(method='ffill'))

def gcdH(park):
    '''
    function to get capacity time series for a wind park
    '''
    cap = windparks[parks==park].Capacity.values
    com = windparks[parks==park].commissioning.astype(np.datetime64).values
    cdf = get_cap_df(cap,com).tz_localize('UTC')
    cdf.name = park
    return(cdf)

def analyse_NZparkh(park):
    '''
    analyse hourly wind power generation for one park
    '''
    comph = pd.DataFrame({'MERRA2':NZm[park],
                         'ERA5':NZe[park],
                         'MERRA2_GWA2':NZmg2[park],
                         'ERA5_GWA2':NZeg2[park],
                         'MERRA2_GWA3':NZmg3[park],
                         'ERA5_GWA3':NZeg3[park]})
    comph['obs'] = comph.index.map(prod_NZh[park])/1000
    # get capacities
    caph = capdfH[park]
    # calculate capacity factors
    cfh = comph.div(caph,axis=0)
    # remove capacity factors > 1
    cfh = cfh.mask(cfh>1).dropna()
    stat_h = pd.DataFrame({'ERA5':stats(cfh.ERA5,cfh.obs,False),
                           'ERA5_GWA2':stats(cfh.ERA5_GWA2,cfh.obs,False),
                           'ERA5_GWA3':stats(cfh.ERA5_GWA3,cfh.obs,False),
                           'MERRA2':stats(cfh.MERRA2,cfh.obs,False),
                           'MERRA2_GWA2':stats(cfh.MERRA2_GWA2,cfh.obs,False),
                           'MERRA2_GWA3':stats(cfh.MERRA2_GWA3,cfh.obs,False),
                           'obs':[np.nan,np.nan,np.nan,cfh.obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_h.columns = ['param','dataset',park]
    return(stat_h.set_index(['param','dataset']).transpose())

def analyse_NZh():
    '''
    analyse hourly wind power generation for NZ
    '''
    mask = (prod_NZh.notna()*capdfH.notna()).replace(0,np.nan)
    comph = pd.DataFrame({'MERRA2':(NZm*mask).sum(axis=1),
                         'ERA5':(NZe*mask).sum(axis=1),
                         'MERRA2_GWA2':(NZmg2*mask).sum(axis=1),
                         'ERA5_GWA2':(NZeg2*mask).sum(axis=1),
                         'MERRA2_GWA3':(NZmg3*mask).sum(axis=1),
                         'ERA5_GWA3':(NZeg3*mask).sum(axis=1)})
    comph['obs'] = comph.index.map((prod_NZh*mask).sum(axis=1))/1000
    # get capacities
    caph = (capdfH*mask).sum(axis=1)
    # calculate capacity factors
    cfh = comph.div(caph,axis=0)
    # remove capacity factors > 1
    cfh = cfh.mask(cfh>1).dropna()
    stat_h = pd.DataFrame({'ERA5':stats(cfh.ERA5,cfh.obs,False),
                           'ERA5_GWA2':stats(cfh.ERA5_GWA2,cfh.obs,False),
                           'ERA5_GWA3':stats(cfh.ERA5_GWA3,cfh.obs,False),
                           'MERRA2':stats(cfh.MERRA2,cfh.obs,False),
                           'MERRA2_GWA2':stats(cfh.MERRA2_GWA2,cfh.obs,False),
                           'MERRA2_GWA3':stats(cfh.MERRA2_GWA3,cfh.obs,False),
                           'obs':[np.nan,np.nan,np.nan,cfh.obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_h.columns = ['param','dataset','NZ']
    return(stat_h.set_index(['param','dataset']).transpose())

def analyse_NZparkd(park):
    '''
    analyse daily wind power generation for one park
    '''
    mask = (prod_NZh[park].notna()*capdfH[park].notna()).replace(0,np.nan)
    comph = pd.DataFrame({'MERRA2':NZm[park]*mask,
                         'ERA5':NZe[park]*mask,
                         'MERRA2_GWA2':NZmg2[park]*mask,
                         'ERA5_GWA2':NZeg2[park]*mask,
                         'MERRA2_GWA3':NZmg3[park]*mask,
                         'ERA5_GWA3':NZeg3[park]*mask})
    comph['obs'] = comph.index.map(prod_NZh[park]*mask)/1000
    # get capacities and mask
    caph = capdfH[park]*mask
    # aggregate daily
    capd = caph.resample('D').sum()
    compd = comph.resample('D').sum()
    # calculate capacity factors
    cfd = compd.div(capd,axis=0)
    # remove capacity factors > 1
    cfd = cfd.mask(cfd>1).dropna()
    stat_d = pd.DataFrame({'ERA5':stats(cfd.ERA5,cfd.obs,False),
                           'ERA5_GWA2':stats(cfd.ERA5_GWA2,cfd.obs,False),
                           'ERA5_GWA3':stats(cfd.ERA5_GWA3,cfd.obs,False),
                           'MERRA2':stats(cfd.MERRA2,cfd.obs,False),
                           'MERRA2_GWA2':stats(cfd.MERRA2_GWA2,cfd.obs,False),
                           'MERRA2_GWA3':stats(cfd.MERRA2_GWA3,cfd.obs,False),
                           'obs':[np.nan,np.nan,np.nan,cfd.obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_d.columns = ['param','dataset',park]
    return(stat_d.set_index(['param','dataset']).transpose())

def analyse_NZd():
    '''
    analyse daily wind power generation for NZ
    '''
    # mask for masking simulated data and capacities
    # (to only use timespans where also observed data are available)
    mask = (prod_NZh.notna()*capdfH.notna()).replace(0,np.nan)
    # mask and aggregate simulated data
    comph = pd.DataFrame({'MERRA2':(NZm*mask).sum(axis=1),
                         'ERA5':(NZe*mask).sum(axis=1),
                         'MERRA2_GWA2':(NZmg2*mask).sum(axis=1),
                         'ERA5_GWA2':(NZeg2*mask).sum(axis=1),
                         'MERRA2_GWA3':(NZmg3*mask).sum(axis=1),
                         'ERA5_GWA3':(NZeg3*mask).sum(axis=1)})
    comph['obs'] = comph.index.map(prod_NZh.sum(axis=1))/1000
    # mask and aggregate capacities
    caph = (capdfH*mask).sum(axis=1)
    # aggregate daily
    compd = comph.resample('D').sum()
    capd = caph.resample('D').sum()
    # calculate capacity factors
    cfdu = compd.div(capd,axis=0).dropna()
    cfd = cfdu.mask(cfdu>1).dropna()
    stat_d = pd.DataFrame({'ERA5':stats(cfd.ERA5,cfd.obs,False),
                           'ERA5_GWA2':stats(cfd.ERA5_GWA2,cfd.obs,False),
                           'ERA5_GWA3':stats(cfd.ERA5_GWA3,cfd.obs,False),
                           'MERRA2':stats(cfd.MERRA2,cfd.obs,False),
                           'MERRA2_GWA2':stats(cfd.MERRA2_GWA2,cfd.obs,False),
                           'MERRA2_GWA3':stats(cfd.MERRA2_GWA3,cfd.obs,False),
                           'obs':[np.nan,np.nan,np.nan,cfd.obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_d.columns = ['param','dataset','NZ']
    return(stat_d.set_index(['param','dataset']).transpose())

def analyse_NZparkm(park):
    '''
    analyse monthly wind power generation for one park
    '''
    # mask for masking simulated data and capacities
    # (to only use timespans where also observed data are available)
    mask = (prod_NZh[park].notna()*capdfH[park].notna()).replace(0,np.nan)
    comph = pd.DataFrame({'MERRA2':NZm[park]*mask,
                         'ERA5':NZe[park]*mask,
                         'MERRA2_GWA2':NZmg2[park]*mask,
                         'ERA5_GWA2':NZeg2[park]*mask,
                         'MERRA2_GWA3':NZmg3[park]*mask,
                         'ERA5_GWA3':NZeg3[park]*mask})
    comph['obs'] = comph.index.map(prod_NZh[park]*mask)/1000
    # get capacities and mask
    caph = capdfH[park]*mask
    # aggregate monthly
    capm = caph.resample('M').sum()
    compm = comph.resample('M').sum()
    # calculate capacity factors
    cfm = compm.div(capm,axis=0)
    # remove capacity factors > 1
    cfm = cfm.mask(cfm>1).dropna()
    stat_m = pd.DataFrame({'ERA5':stats(cfm.ERA5,cfm.obs,False),
                           'ERA5_GWA2':stats(cfm.ERA5_GWA2,cfm.obs,False),
                           'ERA5_GWA3':stats(cfm.ERA5_GWA3,cfm.obs,False),
                           'MERRA2':stats(cfm.MERRA2,cfm.obs,False),
                           'MERRA2_GWA2':stats(cfm.MERRA2_GWA2,cfm.obs,False),
                           'MERRA2_GWA3':stats(cfm.MERRA2_GWA3,cfm.obs,False),
                           'obs':[np.nan,np.nan,np.nan,cfm.obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_m.columns = ['param','dataset',park]
    return(stat_m.set_index(['param','dataset']).transpose())

def analyse_NZm():
    '''
    analyse monthly wind power generation for NZ
    '''
    # mask for masking simulated data and capacities
    # (to only use timespans where also observed data are available)
    mask = (prod_NZh.notna()*capdfH.notna()).replace(0,np.nan)
    # mask and aggregate simulated data
    comph = pd.DataFrame({'MERRA2':(NZm*mask).sum(axis=1),
                         'ERA5':(NZe*mask).sum(axis=1),
                         'MERRA2_GWA2':(NZmg2*mask).sum(axis=1),
                         'ERA5_GWA2':(NZeg2*mask).sum(axis=1),
                         'MERRA2_GWA3':(NZmg3*mask).sum(axis=1),
                         'ERA5_GWA3':(NZeg3*mask).sum(axis=1)})
    comph['obs'] = comph.index.map(prod_NZh.sum(axis=1))/1000
    # mask and aggregate capacities
    caph = (capdfH*mask).sum(axis=1)
    # aggregate monthly
    compm = comph.resample('M').sum()
    capm = caph.resample('M').sum()
    # calculate capacity factors
    cfmu = compm.div(capm,axis=0).dropna()
    cfm = cfmu.mask(cfmu>1).dropna()
    stat_m = pd.DataFrame({'ERA5':stats(cfm.ERA5,cfm.obs,False),
                           'ERA5_GWA2':stats(cfm.ERA5_GWA2,cfm.obs,False),
                           'ERA5_GWA3':stats(cfm.ERA5_GWA3,cfm.obs,False),
                           'MERRA2':stats(cfm.MERRA2,cfm.obs,False),
                           'MERRA2_GWA2':stats(cfm.MERRA2_GWA2,cfm.obs,False),
                           'MERRA2_GWA3':stats(cfm.MERRA2_GWA3,cfm.obs,False),
                           'obs':[np.nan,np.nan,np.nan,cfm.obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_m.columns = ['param','dataset','NZ']
    return(stat_m.set_index(['param','dataset']).transpose())

def tidy_res(results,temp,scale):
    '''
    function for preparing results for merging
     results are results in wide format (one region per line)
     temp is temporal resolution (h,d,m)
     scale is spatial resoltion (park or country)
    '''
    rt = results.transpose().reset_index().melt(id_vars=['param','dataset'], 
                                                value_vars=results.index, var_name='location')
    ts = pd.DataFrame({'temp':temp,'scale':scale},
                      index = range(reduce(mul,results.shape)))
    return(pd.concat([rt,ts],axis=1))


## Analysis
# load windpark data
windparks = pd.read_csv(nz_path + "/windparks_NZ.csv", delimiter=';', parse_dates=['commissioning'])
# dictionary for matching data
d = {'twf1':'twf_1',
     'twf2':'twf_2',
     'twf3':'twf_3',
     'Te Apiti':'te_apiti',
     'Te Rere Hau1':'te_rere_hau',
     'Te Rere Hau2':'te_rere_hau',
     'Te Rere Hau3':'te_rere_hau',
     'Te Uku':'te_uku',
     'White Hill Wind Farm':'white_hill',
     'Westwind Wind Farm':'west_wind'}
# match wind park names
parks = (windparks.ProjectName + windparks.stage.astype(str).replace({'0':''})).replace(d).values

# get capacitiy time series for all parks
pu = np.unique(parks)
capdfH = pd.Series(pu,index = pu).apply(gcdH).transpose()

# load simulated data
NZm = load_results('MERRA2','none',parks)
NZmg2 = load_results('MERRA2','GWA2',parks)
NZmg3 = load_results('MERRA2','GWA3',parks)
NZe = load_results('ERA5','none',parks)
NZeg2 = load_results('ERA5','GWA2',parks)
NZeg3 = load_results('ERA5','GWA3',parks)

# analyse results
# hourly
resNZph = pd.concat(pd.Series(np.unique(parks)).apply(analyse_NZparkh).to_list(),axis=0)
resNZh = analyse_NZh()
# daily
resNZpd = pd.concat(pd.Series(np.unique(parks)).apply(analyse_NZparkd).to_list(),axis=0)
resNZd = analyse_NZd()
# monthly
resNZpm = pd.concat(pd.Series(np.unique(parks)).apply(analyse_NZparkm).to_list(),axis=0)
resNZm = analyse_NZm()

# tidy and merge results
rNZph = tidy_res(resNZph,'h','park')
rNZh = tidy_res(resNZh,'h','country')
rNZpd = tidy_res(resNZpd,'d','park')
rNZd = tidy_res(resNZd,'d','country')
rNZpm = tidy_res(resNZpm,'m','park')
rNZm = tidy_res(resNZm,'m','country')
results = pd.concat([rNZph,rNZh,rNZpd,rNZd,rNZpm,rNZm],axis=0)
results['ds'] = results.dataset.str.split('_').apply(lambda x: x[0])
results['GWA'] = (results.dataset.str.split('_')+['none']).apply(lambda x: x[1])
# save results
results.to_csv(results_path + '/statNZ.csv')