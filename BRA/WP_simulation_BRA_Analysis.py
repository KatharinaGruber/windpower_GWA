import datetime
from dateutil.relativedelta import *
from fuzzywuzzy import fuzz
import glob
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp
import sys
import xarray as xr


from paths_bra import *

sys.path.append('./..')
from refuelplot import *
setup()

from utils import *

gen_path = bra_path + '/generation'


# load generation data
# load usinas hourly
if gen_path + '/hourly/usinas.pkl' not in glob.glob(gen_path + '/hourly/*.pkl'):
    USIh = pd.read_csv(gen_path + '/hourly/Comparativo_Geração_de_Energia_Semana_data_usinas.csv',
                       sep = ';', index_col = 0, parse_dates = True, dayfirst = True).iloc[1:,[6,8]].sort_index()
    # remove missing values
    USIh = USIh.loc[USIh.index.notnull()].dropna()
    USIh.columns = ['usina','prod_GWh']
    # in RIO DO FOGO there is one duplicate hour after one missing hour -> change timestamps of those hours
    idxUSIh = USIh.index.values
    midxUSIh = USIh.reset_index().set_index(['usina','Data Escala de Tempo 1 GE Comp 3']).index
    idxUSIh[midxUSIh.duplicated(keep='last')]  = idxUSIh[midxUSIh.duplicated(keep='first')] - np.timedelta64(1,'h')
    USIh.index = pd.DatetimeIndex(idxUSIh)
    USIhs = USIh.reset_index().set_index(['usina','index']).unstack(level=0).prod_GWh
    USIhs.to_csv(gen_path + '/hourly/usinas.csv')
    USIhs.to_pickle(gen_path + '/hourly/usinas.pkl')
wpUSIhs = pd.read_pickle(gen_path + '/hourly/usinas.pkl')

# load and match aneel and ons windparks
def get_cap_df(cap,comdate):
    com = pd.DataFrame({'capacity': cap}).groupby(comdate).sum()
    cap_cum = com.capacity.cumsum()
    # if only years given for commissioning dates -> gradual capacity increase over year, full capacity at end of year
    if type(cap_cum.index.values[0]) == np.int64:
        cap_cum.index = [np.datetime64(str(int(year))+"-12-31 23:00:00") for year in cap_cum.index.values]
        # create yearly dates at yearends
        drcc = pd.date_range(np.datetime64('2005-12-31 23:00:00'),
                             np.datetime64('2019-12-31 23:00:00'),freq= 'y')
        cap_cum = pd.Series(drcc.map(cap_cum),index = drcc)
        # if first year emtpy: either year before or 0 if nothing before
        if(sum(com.index<2000) > 0):
            cap_cum[0] = com.cumsum()[com.index<2000].max()
        else:
            cap_cum[0] = 0
        # if missing years -> put capacity of year before
        cap_cum = cap_cum.ffill()
    dr = pd.date_range('1/1/2006','31/12/2019 23:00:00',freq = 'h')
    cap_ts = pd.Series(dr.map(cap_cum),index = dr)
    cap_ts[0] = cap_cum[cap_cum.index<=pd.Timestamp('2006-01-01')].max()
    if type(comdate[0]) == np.int64:
        return(cap_ts.interpolate(method='linear'))
    else:
        return(cap_ts.fillna(method='ffill'))

def matchWords(word, statements):
    # function to match a word to different statements
    # output: ratio of matching (0-100) for all provided statements
    results = []
    for s in statements:
        r = fuzz.ratio(word, s)
        results.append(r)
    return results

def match_string(string, array):
    # function for matching casefolded strings
    Slc = string.strip().casefold()
    Alc = [arr.casefold() for arr in array.str.strip().unique()]
    scores = matchWords(Slc, Alc)
    mscore = max(scores)
    strarr = array.unique()[np.where(np.array(scores)==mscore)][0]
    return(string,strarr,mscore)

def match_anl(string):
    # function to match ONS to ANL windparks
    return(match_string(string,ANL2.name))

# load ANEEL and ONS windparks
ONS = pd.read_csv(bra_path + '/ONS_windparks.csv', index_col = 0)
# remove those with CONJUNTO EOLICO - they're there twice and capacities don't match with ANEEL data
ONS = ONS[~ONS.usina.str.contains('CONJUNTO EOLICO')]
# remove some other duplicate windparks
ONS = ONS[[d not in [' CANOA QUEBRADA (E-RV-ACEP)',' PV DO NORDESTE',' SM (SANTA MARIA)',' SÃO BENTO NORTE II'] for d in ONS.usina]]
ANL = pd.read_csv(bra_path + '/turbine_data.csv', index_col = 0)

# characters and strings to replace for better matching
letters = {'Ãµ':'õ',
           'ó':'o',
           'ã':'a',
           'á':'a',
           'â':'a',
           'é':'e',
           'Ã':'A',
           'Á':'A',
           'Â':'A',
           'Ó':'O',
           'É':'E',
           'ú':'u',
           'ô':'o',
           'Ô':'O',
           'ú':'u',
           'Ú':'U',
           'ç':'c',
           'Ç':'C',
           'í':'i',
           'Í':'I',
           'Ê':'E'}
remove = {' 2LER':'',
          ' 2LFA':'',
          ' LFA':'',
          'EOL ':'',
          ' 3LER':'',
          'Usina Eolica ':'',
          'Eólica ':'',
          ' ENERGIAS RENOVAVEIS':'',
#          ' CONJUNTO EOLICO':'',
          '\(E-BV-ACEP\)':'',
          '\(E-RV-ACEP\)':'',
          '\(BELA BISTA\)':'',
          '\(ENERGEN\)':'',
          '\(Antiga Ventos Maranhenses 05\)':'',
          'PARQUE EOLICO ':'',
          ' - N HORIZ':'',
          'ENERGETICA S/A':'',
          '\(ILHEUS\)':'',
          ' EOLOS':'',
          'S\.A\.':''}
replace = {'LAG DO':'LAGOA DO',
           'VENTOS S VICENTE':'VENTOS DE SAO VICENTE',
           'SERRA BABILONIA':'SERRA DA BABILONIA',
           'CORREDOR SENANDES':'CORREDOR DO SENANDES',
           'SAO BENTO NORTE':'SAO BENTO DO NORTE',
           'GAMELEIRAS':'GAMELERIAS',
           'Lagoinha':'Lagoinh',
           'PAPAGAIOS':'PAPAGAIO',
           'VENTOS DE SAO ABRAAO':'VENTOS DO SANTO ABRAAO',
           'VENTOS DO SAO MARIO':'VENTOS DE SAO MARIO',
           'DAGUA':'D AGUA',
           'B VEN':'BONS VENTOS',
           'NOVA BURITI':'BURITI',
           'NOVA CAJUCOCO':'CAJUCOCO',
           'PALMAS':'DE PALMAS',
           'DE PALMARES':'PALMARES',
           'PV DO NORDESTE':'VENTOS DO NORDESTE',
           'Aura Lagoa do Barro':'Lagoa do Barro',
           'AURA LAGOA DO BARRO':'LAGOA DO BARRO',
           'LAGOA BARRO':'LAGOA DO BARRO',
           'GRAVATA':'GRAVATA FRUITRADE',
           'FAZENDA DO ROSARIO':'FAZENDA ROSARIO',
           'Parque Eolico do Horizonte':'Ventos de Horizonte',
           'S BENTO':'SAO BENTO',
           'SANTO ANTONIO (BTG PACTUAL)':'SANTO ANTONIO DE PADUA',
           'SM \(SANTA MARIA\)':'SANTA MARIA',
           'SAO JORGE CE':'SAO JORGE',
           'VENT DA ST ESPERANCA':'VENTOS DA SANTA ESPERANCA',
           'VENTOS DA STA DULCE':'VENTOS DA SANTA DULCE',
           'ESPERANCA NORDESTE':'ESPERANCA DO NORDESTE',
           'Eolica Delta':'Delta',
           'Eolica Serra das Vacas':'Serra das Vacas',
           'Ventos de Santo Augusto':'Santo Augusto',
           'Ventos do Sao Gabriel':'Sao Gabriel',
           'GE Maria Helena':'Maria Helena'}
numbers = {'10':'X',
           '11':'XI',
           '12':'XII',
           '13':'XIII',
           '14':'XIV',
           '15':'XV',
           '17':'XVII',
           '19':'XIX',
           '21':'XXI',
           '23':'XXIII',
           '24':'XXIV',
           '25':'XXV',
           '26':'XXVI',
           '27':'XXVII',
           '28':'XXVIII',
           '29':'XXIX',
           '31':'XXXI',
           '34':'XXXIV',
           '35':'XXXV',
           '36':'XXXVI',
           '01':'I',
           '02':'II',
           '03':'III',
           '04':'IV',
           '05':'V',
           '06':'VI',
           '07':'VII',
           '08':'VIII',
           '09':'IX',
           '1':'I',
           '2':'II',
           '3':'III',
           '4':'IV',
           '5':'V',
           '6':'VI',
           '7':'VII',
           '8':'VIII',
           '9':'IX'}

# replace characters
ONS2 = ONS.copy(deep=True)
ANL2 = ANL.copy(deep=True)
for i in letters:
    ONS2.usina = ONS2.usina.str.replace(i,letters.get(i))
    ANL2.name = ANL2.name.str.replace(i,letters.get(i))
for i in replace:
    ONS2.usina = ONS2.usina.str.replace(i,replace.get(i))
    ANL2.name = ANL2.name.str.replace(i,replace.get(i))
for i in remove:
    ONS2.usina = ONS2.usina.str.replace(i,remove.get(i))
for i in numbers:
    ONS2.usina = ONS2.usina.str.replace(i,numbers.get(i))
    ANL2.name = ANL2.name.str.replace(i,numbers.get(i))

# match windparks
matches = ONS2.usina.apply(match_anl).apply(pd.Series)
matches.columns = ['ONS_name','ANL_name','score']

ONSd = pd.Series(ONS.usina.values,index=ONS2.usina.values)
ANLd = pd.Series(ANL.name.values,index=ANL2.name.values)
ONSd.columns = ['simpl','orig']
ANLd.columns = ['simpl','orig']


# load simulated data
# prepare simulated data as dataframe
if (results_path + '/wpUSI_MER.pkl' not in glob.glob(results_path + '/*.pkl')):
    wpERAxr = xr.open_dataset(results_path + '/windpower_stat_ERA5.nc',chunks={'time':80})
    wpMERxr = xr.open_dataset(results_path + '/windpower_stat_MERRA2.nc',chunks={'time':80})
    wpERAgxr = xr.open_mfdataset(results_path +'/windpower_??_ERA5_GWA.nc',chunks={'time':80})
    wpMERgxr = xr.open_mfdataset(results_path +'/windpower_??_MERRA2_GWA.nc',chunks={'time':80})
    
    turb_mer = pd.read_csv(bra_path + '/turbine_data_mer.csv',index_col=0)
    turb_era = pd.read_csv(bra_path + '/turbine_data_era.csv',index_col=0)
    turb_merg = pd.read_csv(bra_path + '/turbine_data_mer_gwa3.csv',index_col=0)
    turb_erag = pd.read_csv(bra_path + '/turbine_data_era_gwa3.csv',index_col=0)
    
    lbl = pd.read_csv(bra_path+ '/labels_turbine_data_gwa3.csv',index_col=0)
    
    wpMERdf = wpMERxr.to_dataframe().unstack().wp
    wpERAdf = wpERAxr.to_dataframe().unstack().wp
    wpMERgdf = wpMERgxr.assign_coords(location=range(len(wpMERgxr.location.values))).to_dataframe().unstack().wp
    wpERAgdf = wpERAgxr.assign_coords(location=range(len(wpERAgxr.location.values))).to_dataframe().unstack().wp
    # some locations have more than one park, get shares of parks
    sharesMER = ANL.cap.groupby([lbl.lbl_mer.values,ANL.name.values]).sum()/ANL.cap.groupby([lbl.lbl_mer.values,ANL.name.values]).sum().index.get_level_values(0).map(ANL.cap.groupby(lbl.lbl_mer.values).sum())
    sharesERA = ANL.cap.groupby([lbl.lbl_era.values,ANL.name.values]).sum()/ANL.cap.groupby([lbl.lbl_era.values,ANL.name.values]).sum().index.get_level_values(0).map(ANL.cap.groupby(lbl.lbl_era.values).sum())
    sharesMERg = ANL.cap.groupby([lbl.lbl_mer_gwa.values,ANL.name.values]).sum()/ANL.cap.groupby([lbl.lbl_mer_gwa.values,ANL.name.values]).sum().index.get_level_values(0).map(ANL.cap.groupby(lbl.lbl_mer_gwa.values).sum())
    sharesERAg = ANL.cap.groupby([lbl.lbl_era_gwa.values,ANL.name.values]).sum()/ANL.cap.groupby([lbl.lbl_era_gwa.values,ANL.name.values]).sum().index.get_level_values(0).map(ANL.cap.groupby(lbl.lbl_era_gwa.values).sum())
    # get generation per park
    wpMER = wpMERdf.loc[sharesMER.index.codes[0].values()].mul(sharesMER.values,axis=0).groupby(sharesMER.index.get_level_values(1).values).sum().transpose()
    wpERA = wpERAdf.loc[sharesERA.index.codes[0].values()].mul(sharesERA.values,axis=0).groupby(sharesERA.index.get_level_values(1).values).sum().transpose()
    wpMERg = wpMERgdf.loc[sharesMERg.index.codes[0].values()].mul(sharesMERg.values,axis=0).groupby(sharesMERg.index.get_level_values(1).values).sum().transpose()
    wpERAg = wpERAgdf.loc[sharesERAg.index.codes[0].values()].mul(sharesERAg.values,axis=0).groupby(sharesERAg.index.get_level_values(1).values).sum().transpose()
    # adapt index of MERRA data in 2019 (substract half an hour)
    wpMER.index = wpMER.index[wpMER.index<'2018-12'].append(wpMER.index[wpMER.index>='2018-12'] - np.timedelta64(30,'m'))
    wpMERg.index = wpMER.index[wpMERg.index<'2018-12'].append(wpMERg.index[wpMERg.index>='2018-12'] - np.timedelta64(30,'m'))
    # set time zones
    wpMER = wpMER.tz_localize('UTC').tz_convert('America/Bahia')
    wpERA = wpERA.tz_localize('UTC').tz_convert('America/Bahia')
    wpMERg = wpMERg.tz_localize('UTC').tz_convert('America/Bahia')
    wpERAg = wpERAg.tz_localize('UTC').tz_convert('America/Bahia')
    
    wpMER.to_pickle(results_path + '/wpUSI_MER.pkl')
    wpERA.to_pickle(results_path + '/wpUSI_ERA.pkl')
    wpMERg.to_pickle(results_path + '/wpUSI_MERgwa.pkl')
    wpERAg.to_pickle(results_path + '/wpUSI_ERAgwa.pkl')
else:
    wpMER = pd.read_pickle(results_path + '/wpUSI_MER.pkl')
    wpERA = pd.read_pickle(results_path + '/wpUSI_ERA.pkl')
    wpMERg = pd.read_pickle(results_path + '/wpUSI_MERgwa.pkl')
    wpERAg = pd.read_pickle(results_path + '/wpUSI_ERAgwa.pkl')
	
	
# data cleaning
# 0. remove leading and trailing 0s in observed data (replace by nans)
wpUSIhs[wpUSIhs.fillna(0).cumsum(axis=0)==0] = np.nan # remove leading 0s
wpUSIhs[wpUSIhs[::-1].fillna(0).cumsum(axis=0)[::-1]==0] = np.nan # remove trailing 0s
# 1. get  matching power generation timeseries 
matches2 = pd.DataFrame({'ANL_name':matches.ANL_name.map(ANLd.drop_duplicates()),
                         'ONS_name':matches.ONS_name.map(ONSd),
                         'score':matches.score}) # put ANEEL and ONS together into dataframe
matches2H = matches2.copy(deep=True) # get matching hourly windparks
matches2H = matches2H[[usi in wpUSIhs.columns.values for usi in matches2.ONS_name]]
matches2.to_pickle(bra_path + '/matches2.pkl') # save matches
matches2H.to_pickle(bra_path + '/matches2H.pkl')
# 2. remove constant timeseries from observed data
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
    #return(len(sd))
    rmdf = pd.Series(0,index=wpt.index)
    for i in range(len(sd)):
        rmdf.iloc[sd[i]:ed[i]] = 1
    return(wpt.where(rmdf==0))
wpUSIh = wpUSIhs.apply(rm_constTS,axis=0).tz_localize('America/Bahia',ambiguous=False) # part of data removed in 150 of 174 windparks
# 3. separate short time series
# get short time series (less than 2 years)
USIhs = wpUSIh[wpUSIh.columns[wpUSIh.notna().sum(axis=0)<8760*2].values]
matches2Hs = matches2H[[usi in USIhs.columns.values for usi in matches2H.ONS_name.values]]
# get long time series (at least 2 years)
USIhl = wpUSIh.drop(wpUSIh.columns[wpUSIh.notna().sum(axis=0)<8760*2].values,axis=1)
matches2Hl = matches2H[[usi in USIhl.columns.values for usi in matches2H.ONS_name.values]]
# 4. find CFs > 1
def gcdH(park):
    cap = ANL[ANL.name==park].cap.values
    com = ANL[ANL.name==park].commissioning.astype(np.datetime64).values
    return(get_cap_df(cap,com).tz_localize('UTC').tz_convert('America/Bahia'))
# prepare/load capacities
if gen_path+'/usi_capH_filterONS.pkl' not in glob.glob(gen_path + '/*'):
    t1 = datetime.datetime.now()
    capdf = pd.Series(ANL.name.unique()).apply(gcdH)
    t2 = datetime.datetime.now()
    print(t2-t1)
    capdf.index = ANL.name.unique()
    capdfH = capdf.transpose()
    capdfH.to_pickle(gen_path+'/usi_capH_filterONS.pkl')
capdfH = pd.read_pickle(gen_path+'/usi_capH_filterONS.pkl')
# get capacities of matching time series
capUSIhs = capdfH[matches2Hs.ANL_name]
capUSIhl = capdfH[matches2Hl.ANL_name]
# adapt names
capUSIhs.columns = matches2Hs.ONS_name.values
capUSIhl.columns = matches2Hl.ONS_name.values
# calculate capacity factors
cfUSIhs = (USIhs[matches2Hs.ONS_name]/(capUSIhs/10**6))
cfUSIhl = (USIhl[matches2Hl.ONS_name]/(capUSIhl/10**6))
# 12 parks where there are more than 10 days of CFs > 1 in hourly data (11 with matching score 100)
# remove those bad parks and for other with just a few CFs > 1 replace those CFs with nans (happens in statistics calculation function)
USIhlc = USIhl.drop(cfUSIhl.columns[((cfUSIhl>1).sum(axis=0)>10*24)].values,axis=1)
matches2Hlc = matches2Hl.set_index('ONS_name').drop(cfUSIhl.columns[((cfUSIhl>1).sum(axis=0)>10*24)].values).reset_index() # remove also from matching data


# calculate statistics
# Usinas
def analyseUSIh(parks):
    compUSIh= pd.DataFrame({'MERRA2':wpMER[parks.ANL_name],
                            'ERA5':wpERA[parks.ANL_name],
                            'MERRA2_GWA':wpMERg[parks.ANL_name],
                            'ERA5_GWA':wpERAg[parks.ANL_name],
                            'wp_obs':wpUSIh[parks.ONS_name]*10**6})
    # get capacities
    capUSIh = capdfH[parks.ANL_name]
    # calculate capacity factors
    cf_USIh = compUSIh.div(capUSIh,axis=0)
    # remove capacity factors > 1
    cf_USIh = cf_USIh.mask(cf_USIh>1).dropna()
    stat_h = pd.DataFrame({'ERA5':stats(cf_USIh.ERA5,cf_USIh.wp_obs,False),
                           'ERA5_GWA':stats(cf_USIh.ERA5_GWA,cf_USIh.wp_obs,False),
                           'MERRA2':stats(cf_USIh.MERRA2,cf_USIh.wp_obs,False),
                           'MERRA2_GWA':stats(cf_USIh.MERRA2_GWA,cf_USIh.wp_obs,False),
                           'obs':[np.nan,np.nan,np.nan,cf_USIh.wp_obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_h.columns = ['param','dataset',parks.ANL_name]
    return(stat_h.set_index(['param','dataset']).transpose())
def analyseUSId(parks):
    compUSIh= pd.DataFrame({'MERRA2':wpMER[parks.ANL_name],
                            'ERA5':wpERA[parks.ANL_name],
                            'MERRA2_GWA':wpMERg[parks.ANL_name],
                            'ERA5_GWA':wpERAg[parks.ANL_name],
                            'wp_obs':wpUSIh[parks.ONS_name]*10**6})
    # get capacities and mask
    capUSIh = capdfH[parks.ANL_name].where(compUSIh.wp_obs.notna(),np.nan)
    # mask
    compUSIhm = compUSIh.where(compUSIh.wp_obs.notna(),np.nan,axis=1)
    # aggregate daily
    capUSId = capUSIh.resample('D').sum()
    compUSId = compUSIhm.resample('D').sum()
    # calculate capacity factors
    cf_USId = compUSId.div(capUSId,axis=0)
    # remove capacity factors > 1
    cf_USId = cf_USId.mask(cf_USId>1).dropna()
    stat_d = pd.DataFrame({'ERA5':stats(cf_USId.ERA5,cf_USId.wp_obs,False),
                           'ERA5_GWA':stats(cf_USId.ERA5_GWA,cf_USId.wp_obs,False),
                           'MERRA2':stats(cf_USId.MERRA2,cf_USId.wp_obs,False),
                           'MERRA2_GWA':stats(cf_USId.MERRA2_GWA,cf_USId.wp_obs,False),
                           'obs':[np.nan,np.nan,np.nan,cf_USId.wp_obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_d.columns = ['param','dataset',parks.ANL_name]
    return(stat_d.set_index(['param','dataset']).transpose())
def analyseUSIm(parks):
    compUSIh= pd.DataFrame({'MERRA2':wpMER[parks.ANL_name],
                            'ERA5':wpERA[parks.ANL_name],
                            'MERRA2_GWA':wpMERg[parks.ANL_name],
                            'ERA5_GWA':wpERAg[parks.ANL_name],
                            'wp_obs':wpUSIh[parks.ONS_name]*10**6})
    # get capacities and mask
    capUSIh = capdfH[parks.ANL_name].where(compUSIh.wp_obs.notna(),np.nan)
    # mask
    compUSIhm = compUSIh.where(compUSIh.wp_obs.notna(),np.nan,axis=1)
    # aggregate monthly
    capUSIm = capUSIh.resample('M').sum()
    compUSIm = compUSIhm.resample('M').sum()
    # calculate capacity factors
    cf_USIm = compUSIm.div(capUSIm,axis=0)
    # remove capacity factors > 1
    cf_USIm = cf_USIm.mask(cf_USIm>1).dropna()
    stat_m = pd.DataFrame({'ERA5':stats(cf_USIm.ERA5,cf_USIm.wp_obs,False),
                           'ERA5_GWA':stats(cf_USIm.ERA5_GWA,cf_USIm.wp_obs,False),
                           'MERRA2':stats(cf_USIm.MERRA2,cf_USIm.wp_obs,False),
                           'MERRA2_GWA':stats(cf_USIm.MERRA2_GWA,cf_USIm.wp_obs,False),
                           'obs':[np.nan,np.nan,np.nan,cf_USIm.wp_obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_m.columns = ['param','dataset',parks.ANL_name]
    return(stat_m.set_index(['param','dataset']).transpose())
stats_USIh = pd.concat(matches2Hlc[matches2Hlc.score==100].apply(analyseUSIh,axis=1).tolist(),axis=0).transpose().dropna(axis=1)
stats_USId = pd.concat(matches2Hlc[matches2Hlc.score==100].apply(analyseUSId,axis=1).tolist(),axis=0).transpose().dropna(axis=1)
stats_USIm = pd.concat(matches2Hlc[matches2Hlc.score==100].apply(analyseUSIm,axis=1).tolist(),axis=0).transpose().dropna(axis=1)
stats_USIh.to_csv(results_path + '/stats_USIh.csv')
stats_USId.to_csv(results_path + '/stats_USId.csv')
stats_USIm.to_csv(results_path + '/stats_USIm.csv')
# States
# insert states in matches dataframes
matches2Hs['state'] = matches2Hs.ANL_name.map(ANL.groupby('name').state.first()).values
matches2Hlc['state'] = matches2Hlc.ANL_name.map(ANL.groupby('name').state.first()).values
matches2Hlc.to_pickle(bra_path + '/matches2Hlc.pkl')
def analyseESTh(usinas):
    # remove leading and trailing 0s in observed data
    wpobs = wpUSIh[usinas.ONS_name.values].sum(axis=1).copy(deep=True)
    wpobs[wpobs.cumsum(axis=0)==0] = np.nan
    wpobs[wpobs[::-1].cumsum(axis=0)[::-1]==0] = np.nan
    # mask for masking simulated data and capacities
    # (to only use timespans where also observed data are available)
    mask = wpUSIh[usinas.ONS_name.values].notna()
    mask.columns = usinas.ANL_name.values
    # mask and aggregate simulated data
    wpMER_ESTh = (wpMER[usinas.ANL_name.values]*mask).sum(axis=1)
    wpERA_ESTh = (wpERA[usinas.ANL_name.values]*mask).sum(axis=1)
    wpMERg_ESTh = (wpMERg[usinas.ANL_name.values]*mask).sum(axis=1)
    wpERAg_ESTh = (wpERAg[usinas.ANL_name.values]*mask).sum(axis=1)
    compESTh= pd.DataFrame({'MERRA2':wpMER_ESTh,
                            'ERA5':wpERA_ESTh,
                            'MERRA2_GWA':wpMERg_ESTh,
                            'ERA5_GWA':wpERAg_ESTh,
                            'wp_obs':wpobs*10**6})
    # mask and aggregate capacities
    capusish = capdfH[usinas.ANL_name.values]*mask
    capESTh = capusish.sum(axis=1)
    # calculate capacity factors
    cf_ESThu = compESTh.div(capESTh,axis=0).dropna()
    cf_ESTh = cf_ESThu.mask(cf_ESThu>1).dropna()
    stat_h = pd.DataFrame({'ERA5':stats(cf_ESTh.ERA5,cf_ESTh.wp_obs,False),
                           'ERA5_GWA':stats(cf_ESTh.ERA5_GWA,cf_ESTh.wp_obs,False),
                           'MERRA2':stats(cf_ESTh.MERRA2,cf_ESTh.wp_obs,False),
                           'MERRA2_GWA':stats(cf_ESTh.MERRA2_GWA,cf_ESTh.wp_obs,False),
                           'obs':[np.nan,np.nan,np.nan,cf_ESTh.wp_obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_h.columns = ['param','dataset',usinas.state.values[0]]
    return(stat_h.set_index(['param','dataset']).transpose())
def analyseESTd(usinas):
    # remove leading and trailing 0s in observed data
    wpobs = wpUSIh[usinas.ONS_name.values].sum(axis=1).copy(deep=True)
    wpobs[wpobs.cumsum(axis=0)==0] = np.nan
    wpobs[wpobs[::-1].cumsum(axis=0)[::-1]==0] = np.nan
    # mask for masking simulated data and capacities
    # (to only use timespans where also observed data are available)
    mask = wpUSIh[usinas.ONS_name.values].notna()
    mask.columns = usinas.ANL_name.values
    # mask and aggregate simulated data
    wpMER_ESTh = (wpMER[usinas.ANL_name.values]*mask).sum(axis=1)
    wpERA_ESTh = (wpERA[usinas.ANL_name.values]*mask).sum(axis=1)
    wpMERg_ESTh = (wpMERg[usinas.ANL_name.values]*mask).sum(axis=1)
    wpERAg_ESTh = (wpERAg[usinas.ANL_name.values]*mask).sum(axis=1)
    compESTh= pd.DataFrame({'MERRA2':wpMER_ESTh,
                            'ERA5':wpERA_ESTh,
                            'MERRA2_GWA':wpMERg_ESTh,
                            'ERA5_GWA':wpERAg_ESTh,
                            'wp_obs':wpobs*10**6})
    # mask and aggregate capacities
    capusish = capdfH[usinas.ANL_name.values]*mask
    capESTh = capusish.sum(axis=1)
    # aggregate daily
    compESTd = compESTh.resample('D').sum()
    capESTd = capESTh.resample('D').sum()
    # calculate capacity factors
    cf_ESTdu = compESTd.div(capESTd,axis=0).dropna()
    cf_ESTd = cf_ESTdu.mask(cf_ESTdu>1).dropna()
    stat_d = pd.DataFrame({'ERA5':stats(cf_ESTd.ERA5,cf_ESTd.wp_obs,False),
                           'ERA5_GWA':stats(cf_ESTd.ERA5_GWA,cf_ESTd.wp_obs,False),
                           'MERRA2':stats(cf_ESTd.MERRA2,cf_ESTd.wp_obs,False),
                           'MERRA2_GWA':stats(cf_ESTd.MERRA2_GWA,cf_ESTd.wp_obs,False),
                           'obs':[np.nan,np.nan,np.nan,cf_ESTd.wp_obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_d.columns = ['param','dataset',usinas.state.values[0]]
    return(stat_d.set_index(['param','dataset']).transpose())
def analyseESTm(usinas):
    # remove leading and trailing 0s in observed data
    wpobs = wpUSIh[usinas.ONS_name.values].sum(axis=1).copy(deep=True)
    wpobs[wpobs.cumsum(axis=0)==0] = np.nan
    wpobs[wpobs[::-1].cumsum(axis=0)[::-1]==0] = np.nan
    # mask for masking simulated data and capacities
    # (to only use timespans where also observed data are available)
    mask = wpUSIh[usinas.ONS_name.values].notna()
    mask.columns = usinas.ANL_name.values
    # mask and aggregate simulated data
    wpMER_ESTh = (wpMER[usinas.ANL_name.values]*mask).sum(axis=1)
    wpERA_ESTh = (wpERA[usinas.ANL_name.values]*mask).sum(axis=1)
    wpMERg_ESTh = (wpMERg[usinas.ANL_name.values]*mask).sum(axis=1)
    wpERAg_ESTh = (wpERAg[usinas.ANL_name.values]*mask).sum(axis=1)
    compESTh= pd.DataFrame({'MERRA2':wpMER_ESTh,
                            'ERA5':wpERA_ESTh,
                            'MERRA2_GWA':wpMERg_ESTh,
                            'ERA5_GWA':wpERAg_ESTh,
                            'wp_obs':wpobs*10**6})
    # mask and aggregate capacities
    capusish = capdfH[usinas.ANL_name.values]*mask
    capESTh = capusish.sum(axis=1)
    # aggregate monthly
    compESTm = compESTh.resample('M').sum()
    capESTm = capESTh.resample('M').sum()
    # calculate capacity factors
    cf_ESTmu = compESTm.div(capESTm,axis=0).dropna()
    cf_ESTm = cf_ESTmu.mask(cf_ESTmu>1).dropna()
    stat_m = pd.DataFrame({'ERA5':stats(cf_ESTm.ERA5,cf_ESTm.wp_obs,False),
                           'ERA5_GWA':stats(cf_ESTm.ERA5_GWA,cf_ESTm.wp_obs,False),
                           'MERRA2':stats(cf_ESTm.MERRA2,cf_ESTm.wp_obs,False),
                           'MERRA2_GWA':stats(cf_ESTm.MERRA2_GWA,cf_ESTm.wp_obs,False),
                           'obs':[np.nan,np.nan,np.nan,cf_ESTm.wp_obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_m.columns = ['param','dataset',usinas.state.values[0]]
    return(stat_m.set_index(['param','dataset']).transpose())
stats_ESTh = matches2Hlc[matches2Hlc.score==100].groupby('state').apply(analyseESTh).transpose()
stats_ESTh.columns = stats_ESTh.columns.get_level_values(0).values
stats_ESTd = matches2Hlc[matches2Hlc.score==100].groupby('state').apply(analyseESTd).transpose()
stats_ESTd.columns = stats_ESTd.columns.get_level_values(0).values
stats_ESTm = matches2Hlc[matches2Hlc.score==100].groupby('state').apply(analyseESTm).transpose()
stats_ESTm.columns = stats_ESTm.columns.get_level_values(0).values
stats_ESTh.to_csv(results_path + '/stats_ESTh.csv')
stats_ESTd.to_csv(results_path + '/stats_ESTd.csv')
stats_ESTm.to_csv(results_path + '/stats_ESTm.csv')
# Brazil
def analyseBRAh(usinas):
    # remove leading and trailing 0s in observed data
    wpobs = wpUSIh[usinas.ONS_name.values].sum(axis=1).copy(deep=True)
    wpobs[wpobs.cumsum(axis=0)==0] = np.nan
    wpobs[wpobs[::-1].cumsum(axis=0)[::-1]==0] = np.nan
    # mask for masking simulated data and capacities
    # (to only use timespans where also observed data are available)
    mask = wpUSIh[usinas.ONS_name.values].notna()
    mask.columns = usinas.ANL_name.values
    # mask and aggregate simulated data
    wpMER_BRAh = (wpMER[usinas.ANL_name.values]*mask).sum(axis=1)
    wpERA_BRAh = (wpERA[usinas.ANL_name.values]*mask).sum(axis=1)
    wpMERg_BRAh = (wpMERg[usinas.ANL_name.values]*mask).sum(axis=1)
    wpERAg_BRAh = (wpERAg[usinas.ANL_name.values]*mask).sum(axis=1)
    compBRAh= pd.DataFrame({'MERRA2':wpMER_BRAh,
                            'ERA5':wpERA_BRAh,
                            'MERRA2_GWA':wpMERg_BRAh,
                            'ERA5_GWA':wpERAg_BRAh,
                            'wp_obs':wpobs*10**6})
    # mask and aggregate capacities
    capusish = capdfH[usinas.ANL_name.values]*mask
    capBRAh = capusish.sum(axis=1)
    # calculate capacity factors
    cf_BRAhu = compBRAh.div(capBRAh,axis=0).dropna()
    cf_BRAh = cf_BRAhu.mask(cf_BRAhu>1).dropna()
    stat_h = pd.DataFrame({'ERA5':stats(cf_BRAh.ERA5,cf_BRAh.wp_obs,False),
                           'ERA5_GWA':stats(cf_BRAh.ERA5_GWA,cf_BRAh.wp_obs,False),
                           'MERRA2':stats(cf_BRAh.MERRA2,cf_BRAh.wp_obs,False),
                           'MERRA2_GWA':stats(cf_BRAh.MERRA2_GWA,cf_BRAh.wp_obs,False),
                           'obs':[np.nan,np.nan,np.nan,cf_BRAh.wp_obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_h.columns = ['param','dataset','BRA']
    return(stat_h.set_index(['param','dataset']).transpose())
def analyseBRAd(usinas):
    # remove leading and trailing 0s in observed data
    wpobs = wpUSIh[usinas.ONS_name.values].sum(axis=1).copy(deep=True)
    wpobs[wpobs.cumsum(axis=0)==0] = np.nan
    wpobs[wpobs[::-1].cumsum(axis=0)[::-1]==0] = np.nan
    # mask for masking simulated data and capacities
    # (to only use timespans where also observed data are available)
    mask = wpUSIh[usinas.ONS_name.values].notna()
    mask.columns = usinas.ANL_name.values
    # mask and aggregate simulated data
    wpMER_BRAh = (wpMER[usinas.ANL_name.values]*mask).sum(axis=1)
    wpERA_BRAh = (wpERA[usinas.ANL_name.values]*mask).sum(axis=1)
    wpMERg_BRAh = (wpMERg[usinas.ANL_name.values]*mask).sum(axis=1)
    wpERAg_BRAh = (wpERAg[usinas.ANL_name.values]*mask).sum(axis=1)
    compBRAh= pd.DataFrame({'MERRA2':wpMER_BRAh,
                            'ERA5':wpERA_BRAh,
                            'MERRA2_GWA':wpMERg_BRAh,
                            'ERA5_GWA':wpERAg_BRAh,
                            'wp_obs':wpobs*10**6})
    # mask and aggregate capacities
    capusish = capdfH[usinas.ANL_name.values]*mask
    capBRAh = capusish.sum(axis=1)
    # aggregate daily
    compBRAd = compBRAh.resample('D').sum()
    capBRAd = capBRAh.resample('D').sum()
    # calculate capacity factors
    cf_BRAdu = compBRAd.div(capBRAd,axis=0).dropna()
    cf_BRAd = cf_BRAdu.mask(cf_BRAdu>1).dropna()
    stat_d = pd.DataFrame({'ERA5':stats(cf_BRAd.ERA5,cf_BRAd.wp_obs,False),
                           'ERA5_GWA':stats(cf_BRAd.ERA5_GWA,cf_BRAd.wp_obs,False),
                           'MERRA2':stats(cf_BRAd.MERRA2,cf_BRAd.wp_obs,False),
                           'MERRA2_GWA':stats(cf_BRAd.MERRA2_GWA,cf_BRAd.wp_obs,False),
                           'obs':[np.nan,np.nan,np.nan,cf_BRAd.wp_obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_d.columns = ['param','dataset','BRA']
    return(stat_d.set_index(['param','dataset']).transpose())

def analyseBRAm(usinas):
    # remove leading and trailing 0s in observed data
    wpobs = wpUSIh[usinas.ONS_name.values].sum(axis=1).copy(deep=True)
    wpobs[wpobs.cumsum(axis=0)==0] = np.nan
    wpobs[wpobs[::-1].cumsum(axis=0)[::-1]==0] = np.nan
    # mask for masking simulated data and capacities
    # (to only use timespans where also observed data are available)
    mask = wpUSIh[usinas.ONS_name.values].notna()
    mask.columns = usinas.ANL_name.values
    # mask and aggregate simulated data
    wpMER_BRAh = (wpMER[usinas.ANL_name.values]*mask).sum(axis=1)
    wpERA_BRAh = (wpERA[usinas.ANL_name.values]*mask).sum(axis=1)
    wpMERg_BRAh = (wpMERg[usinas.ANL_name.values]*mask).sum(axis=1)
    wpERAg_BRAh = (wpERAg[usinas.ANL_name.values]*mask).sum(axis=1)
    compBRAh= pd.DataFrame({'MERRA2':wpMER_BRAh,
                            'ERA5':wpERA_BRAh,
                            'MERRA2_GWA':wpMERg_BRAh,
                            'ERA5_GWA':wpERAg_BRAh,
                            'wp_obs':wpobs*10**6})
    # mask and aggregate capacities
    capusish = capdfH[usinas.ANL_name.values]*mask
    capBRAh = capusish.sum(axis=1)
    # aggregate monthly
    compBRAm = compBRAh.resample('M').sum()
    capBRAm = capBRAh.resample('M').sum()
    # calculate capacity factors
    cf_BRAmu = compBRAm.div(capBRAm,axis=0).dropna()
    cf_BRAm = cf_BRAmu.mask(cf_BRAmu>1).dropna()
    stat_m = pd.DataFrame({'ERA5':stats(cf_BRAm.ERA5,cf_BRAm.wp_obs,False),
                           'ERA5_GWA':stats(cf_BRAm.ERA5_GWA,cf_BRAm.wp_obs,False),
                           'MERRA2':stats(cf_BRAm.MERRA2,cf_BRAm.wp_obs,False),
                           'MERRA2_GWA':stats(cf_BRAm.MERRA2_GWA,cf_BRAm.wp_obs,False),
                           'obs':[np.nan,np.nan,np.nan,cf_BRAm.wp_obs.mean()]},
                          index = ['cor','rmse','mbe','avg']).reset_index().melt(id_vars=['index']).dropna()
    stat_m.columns = ['param','dataset','BRA']
    return(stat_m.set_index(['param','dataset']).transpose())
stats_BRAh = analyseBRAh(matches2Hlc[matches2Hlc.score==100]).transpose()
stats_BRAd = analyseBRAd(matches2Hlc[matches2Hlc.score==100]).transpose()
stats_BRAm = analyseBRAm(matches2Hlc[matches2Hlc.score==100]).transpose()
stats_BRAh.to_csv(results_path + '/stats_BRAh.csv')
stats_BRAd.to_csv(results_path + '/stats_BRAd.csv')
stats_BRAm.to_csv(results_path + '/stats_BRAm.csv')

# merge results and save
sUSIh = stats_USIh.reset_index().melt(id_vars=['param','dataset'], value_vars=stats_USIh.columns, var_name='location')
sUSIh['temp'] = 'h'
sUSIh['scale'] = 'park'
sUSId = stats_USId.reset_index().melt(id_vars=['param','dataset'], value_vars=stats_USId.columns, var_name='location')
sUSId['temp'] = 'd'
sUSId['scale'] = 'park'
sUSIm = stats_USIm.reset_index().melt(id_vars=['param','dataset'], value_vars=stats_USIm.columns, var_name='location')
sUSIm['temp'] = 'm'
sUSIm['scale'] = 'park'
sESTh = stats_ESTh.reset_index().melt(id_vars=['param','dataset'], value_vars=stats_ESTh.columns, var_name='location')
sESTh['temp'] = 'h'
sESTh['scale'] = 'state'
sESTd = stats_ESTd.reset_index().melt(id_vars=['param','dataset'], value_vars=stats_ESTd.columns, var_name='location')
sESTd['temp'] = 'd'
sESTd['scale'] = 'state'
sESTm = stats_ESTm.reset_index().melt(id_vars=['param','dataset'], value_vars=stats_ESTm.columns, var_name='location')
sESTm['temp'] = 'm'
sESTm['scale'] = 'state'
sBRAh = stats_BRAh.reset_index().melt(id_vars=['param','dataset'], value_vars=stats_BRAh.columns, var_name='location')
sBRAh['temp'] = 'h'
sBRAh['scale'] = 'country'
sBRAd = stats_BRAd.reset_index().melt(id_vars=['param','dataset'], value_vars=stats_BRAd.columns, var_name='location')
sBRAd['temp'] = 'd'
sBRAd['scale'] = 'country'
sBRAm = stats_BRAm.reset_index().melt(id_vars=['param','dataset'], value_vars=stats_BRAm.columns, var_name='location')
sBRAm['temp'] = 'm'
sBRAm['scale'] = 'country'
stats = pd.concat([sUSIh, sESTh, sBRAh,
                   sUSId, sESTd, sBRAd,
                   sUSIm, sESTm, sBRAm], axis=0)
stats['GWA'] = 'GWA3'
stats.to_csv(results_path + '/stats_GWA3.csv')
