import numpy as np
import pandas as pd
import xarray as xr
import sys

from paths import results_path

rp_usa = results_path + '/USA'
rp_usa2 = rp_usa + '/results_GWA2'
rp_bra = results_path + '/BRA'
rp_bra2 = rp_bra + '/results_GWA2'

## Load USA data
# load size indicator
nums_usa = pd.read_csv(rp_usa + '/number_grid_points.csv',index_col=0).drop([2,59]) # drop 2,59 - what is it? new england?
nums_usa['scale'] = 'state'
nums_usa.loc[nums_usa.region=='USA','scale'] = 'country'
nums_usa['country'] = 'USA'
# read results
results_USA_GWA2 = pd.read_csv(rp_usa2 + '/stats_GWA2.csv',index_col=0)
results_USA_GWA3 = pd.read_csv(rp_usa + '/stats_GWA3.csv',index_col=0)
results_USA_tidy = pd.concat([results_USA_GWA2,results_USA_GWA3],axis=0).rename({'location':'region'},axis=1)
results_USA_tidy['country'] = 'USA'
# filter regions with bad observed time series (filtered by visual inspection)
bad_states = ['CT','MA','IL','RI','VT','OH','NJ','DE','NC']
bad_states = ['CT','MA','IL','RI','VT','OH','NJ','DE','NC','NE','MI','WI','TN','ND','SD','AK']
bad_regions = ['NewEng','ESC','PacNon']
# filter bad regions
results_USA_tidy = results_USA_tidy.set_index(['country','region']).drop(list(zip(['USA']*len(bad_regions),bad_regions)),axis=0).reset_index()
# filter bad states
results_USA_tidy = results_USA_tidy.set_index(['country','region']).drop(list(zip(['USA']*len(bad_states),bad_states)),axis=0).reset_index()

## Load BRA data
# load size indicator
nums_bra = pd.read_csv(rp_bra + '/number_grid_points.csv',index_col=0)
nums_bra['country'] = 'BRA'
# read results
results_BRA_GWA2 = pd.read_csv(rp_bra2 + '/stats_GWA2.csv',index_col=0)
results_BRA_GWA3 = pd.read_csv(rp_bra + '/stats_GWA3.csv',index_col=0)
results_BRA_tidy = pd.concat([results_BRA_GWA2,results_BRA_GWA3],axis=0).rename({'location':'region'},axis=1)
results_BRA_tidy['country'] = 'BRA'


## prepare results and size parameter for analysis
# calculate "system size" as number between 0 and 1 by divinding by largest number of grid cells per dataset (country)
ss_usa = nums_usa.copy()
ss_bra = nums_bra.copy()
ss_usa.loc[ss_usa.dataset=='MERRA2','cor'] = nums_usa.cor[nums_usa.dataset=='MERRA2']/nums_usa.cor[nums_usa.dataset=='MERRA2'].max()
ss_usa.loc[ss_usa.dataset=='ERA5','cor'] = nums_usa.cor[nums_usa.dataset=='ERA5']/nums_usa.cor[nums_usa.dataset=='ERA5'].max()
ss_bra.loc[ss_bra.dataset=='MERRA2','cor'] = nums_bra.cor[nums_bra.dataset=='MERRA2']/nums_bra.cor[nums_bra.dataset=='MERRA2'].max()
ss_bra.loc[ss_bra.dataset=='ERA5','cor'] = nums_bra.cor[nums_bra.dataset=='ERA5']/nums_bra.cor[nums_bra.dataset=='ERA5'].max()
## add size parameter to results dataframe
# extract dataset only without GWA
results_BRA_tidy['ds'] = results_BRA_tidy.dataset
results_BRA_tidy['dataset'] = results_BRA_tidy.dataset.str.split('_').apply(lambda x: x[0])
results_USA_tidy['ds'] = results_USA_tidy.dataset
results_USA_tidy['dataset'] = results_USA_tidy.dataset.str.split('_').apply(lambda x: x[0])
# add size to results and use cols as index for matching
cols = ['dataset','region','scale','temp']
results_BRA_tidy['systemsize'] = results_BRA_tidy.set_index(cols).index.map(ss_bra.set_index(cols).cor)
results_USA_tidy['systemsize'] = results_USA_tidy.set_index(cols[:-1]).index.map(ss_usa.set_index(cols[:-1]).cor)

## Merge results
results = pd.concat([results_USA_tidy,
                     results_BRA_tidy],
                    axis=0, sort=True, ignore_index=True)

## Clean up data
# remove gwa from non corrected results and remove duplicates
results.loc[(results.ds!='ERA5_GWA')&(results.ds!='MERRA2_GWA'),'GWA'] = 'none'
results = results[~results.duplicated(subset=['country','ds','param','region','scale','systemsize','temp','value'])]
# add gwa to dataset
results['ds2'] = results.ds
results['ds2'][results.GWA!='none'] = results.dataset[results.GWA!='none'] + '_' + results.GWA[results.GWA!='none']

results.to_pickle(results_path + '/results_merged.pkl')
                    
