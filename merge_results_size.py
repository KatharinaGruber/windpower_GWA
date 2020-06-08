import numpy as np
import pandas as pd
import xarray as xr
import sys

from paths import results_path

rp_usa = results_path + '/USA'
rp_usa2 = rp_usa + '/results_GWA2'
rp_bra = results_path + '/BRA'
rp_bra2 = rp_bra + '/results_GWA2'
rp_nz = results_path + '/NZ'

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
results_USA_tidy.scale = results_USA_tidy.scale.replace({'subsystem':'state'})
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

## Load NZ data
# load size indicator
nums_nz = pd.read_csv(rp_nz + '/number_grid_points.csv',index_col=0)
nums_nz['country'] = 'NZ'
# read results
results_NZ_tidy = pd.read_csv(rp_nz + '/statNZ.csv',index_col=0)
results_NZ_tidy.loc[results_NZ_tidy.GWA != 'none','ds'] = results_NZ_tidy.dataset[results_NZ_tidy.GWA != 'none'] + '_GWA'



## prepare results and size parameter for analysis
def merge_res_ss(sys_size,results,id_cols):
    '''
    function for merging the system size parameter and the results
    
    parameters:
     sys_size: dataframe with system size indicator (column cor)
     results: results of statistical analysis in tidy format
     id_cols: columns in sys_size and results used to match the data
    '''
    # calculate "system size" as number between 0 and 1 by divinding by largest number of grid cells per dataset (country)
    sys_size.loc[sys_size.dataset=='MERRA2','cor'] = sys_size.cor[sys_size.dataset=='MERRA2']/sys_size.cor[sys_size.dataset=='MERRA2'].max()
    sys_size.loc[sys_size.dataset=='ERA5','cor'] = sys_size.cor[sys_size.dataset=='ERA5']/sys_size.cor[sys_size.dataset=='ERA5'].max()
    # add size parameter to results dataframe
    # extract dataset only without GWA
    results['ds'] = results.dataset
    results['dataset'] = results.dataset.str.split('_').apply(lambda x: x[0])
    # add size to results and use cols as index for matching
    results['systemsize'] = results.set_index(id_cols).index.map(sys_size.set_index(id_cols).cor)
    return results

# merge size parameter and results
results_USAss = merge_res_ss(nums_usa,results_USA_tidy,['dataset','region','scale'])
results_BRAss = merge_res_ss(nums_bra,results_BRA_tidy,['dataset','region','scale','temp'])
results_NZss = merge_res_ss(nums_nz,results_NZ_tidy,['dataset','scale','temp'])

## Merge results
results = pd.concat([results_USAss,
                     results_BRAss,
                     results_NZss],
                    axis=0, sort=True, ignore_index=True)

## Clean up data
# remove gwa from non corrected results and remove duplicates
results.loc[(results.ds!='ERA5_GWA')&(results.ds!='MERRA2_GWA'),'GWA'] = 'none'
results = results[~results.duplicated(subset=['country','ds','param','region','scale','systemsize','temp','value'])]
# add gwa to dataset
results['ds2'] = results.ds
results.loc[results.GWA!='none','ds2'] = results.dataset[results.GWA!='none'] + '_' + results.GWA[results.GWA!='none']

results.to_pickle(results_path + '/results_merged.pkl')
                    
