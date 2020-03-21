import argparse
import datetime
import glob
import math
import numpy as np
import os
import pandas as pd
import rasterio
import seaborn as sns
import xarray as xr

from matplotlib import pyplot as plt

from dask.diagnostics import ProgressBar
ProgressBar().register()

from paths_usa import *

parser = argparse.ArgumentParser(description='Insert optionally GWA')
parser.add_argument('-GWA')
args = parser.parse_args()
if(args.GWA == None):
    GWA = 3
else:
    GWA = args.GWA

if GWA == "2":
    results_pathg = results_path + '/results_GWA2'
else:
    results_pathg = results_path

# sum up states
if len(glob.glob(results_pathg + '/windpower_states_*.nc')) < 4:
    print("summing up states...")

    # MERRA2
    ofile = results_path + "/windpower_states_MERRA2.nc"
    if ofile not in glob.glob(results_path + '/*'):
        turbine_data_mer = pd.read_csv(usa_path+"/turbine_data_mer.csv")
        wp_loc = xr.open_dataset(results_path+"/windpower_stat_MERRA2.nc")
        wp = wp_loc.groupby(xr.DataArray(turbine_data_mer.state,dims='location')).sum('location')
        wp.to_netcdf(ofile)
        del(wp,wp_loc,turbine_data_mer)

    # MERRA2 + GWA
    ofile = results_pathg + "/windpower_states_MERRA2_GWA.nc"
    if ofile not in glob.glob(results_pathg + '/*'):
        turbine_data_mer_gwa = pd.read_csv(usa_path+"/turbine_data_mer_gwa.csv")
        wp_loc = xr.open_mfdataset(results_pathg+"/windpower_??_MERRA2_GWA.nc")
        wp = wp_loc.groupby(xr.DataArray(turbine_data_mer_gwa.state,dims='location')).sum('location')
        wp.to_netcdf(ofile)
        del(wp,wp_loc,turbine_data_mer_gwa)

    # ERA5
    ofile = results_path + "/windpower_states_ERA5.nc"
    if ofile not in glob.glob(results_path + '/*'):
        turbine_data_era = pd.read_csv(usa_path+"/turbine_data_era.csv")
        wp_loc = xr.open_dataset(results_path+"/windpower_stat_ERA5.nc")
        wp = wp_loc.groupby(xr.DataArray(turbine_data_era.state,dims='location')).sum('location')
        wp.to_netcdf(ofile)
        del(wp,wp_loc,turbine_data_era)

    # ERA5 + GWA
    ofile = results_pathg + "/windpower_states_ERA5_GWA.nc"
    if ofile not in glob.glob(results_pathg + '/*'):
        turbine_data_era_gwa = pd.read_csv(usa_path+"/turbine_data_era_gwa.csv")
        wp_loc = xr.open_mfdataset(results_pathg+"/windpower_??_ERA5_GWA.nc", chunks = {'time':38})
        wp = wp_loc.groupby(xr.DataArray(turbine_data_era_gwa.state,dims='location')).sum('location')
        wp.to_netcdf(ofile)
        del(wp,wp_loc,turbine_data_era_gwa)


# sum up BPA
# the shares actually would depend on the turbine capacity too, but in this case,
# turbines with the same specific power (and thus same label) also have the same capacity
if len(glob.glob(results_pathg + "/windpower_BPA_*.nc")) < 4:
    print("summing up BPA...")

    # get BPA windparks
    BPA_parks = pd.read_csv(usa_path+"/BPA_windparks.csv")
    # get windturbine locations/names
    windturbines = pd.read_csv(usa_path+"/uswtdb_v2_3_20200109.csv",delimiter=',')
    # get labels
    labels = pd.read_csv(usa_path + '/labels_turbine_data.csv')
    # get indices of BPA wind parks from wind turbine dataset
    pBPA = pd.DataFrame({'p': [park in BPA_parks.name.values for park in windturbines[windturbines.t_state!='GU'].p_name.values]})
    # remove turbines that are in other states - error in data?
    pBPA[windturbines[windturbines.t_state!='GU'].xlong.values<-125] = False
    pBPA[windturbines[windturbines.t_state!='GU'].xlong.values>-115] = False

    # MERRA2
    ofile = results_path + "/windpower_BPA_MERRA2.nc"
    if ofile not in glob.glob(results_path + "/*"):
        BPA_mer = pBPA.groupby(labels.lbl_mer).mean() # get shares of where BPA wind parks are installed
        wp_loc_mer = xr.open_dataset(results_path+"/windpower_stat_MERRA2.nc") # load results
        wp_BPA_mer = (wp_loc_mer * np.array(BPA_mer.p)).sum('location') # aggregate for BPA
        wp_BPA_mer.to_netcdf(ofile)
        del(BPA_mer,wp_loc_mer,wp_BPA_mer)
    
    # MERRA2 + GWA
    ofile = results_pathg + "/windpower_BPA_MERRA2_GWA.nc"
    if ofile not in glob.glob(results_pathg + "/*"):
        BPA_mer_gwa = pBPA.groupby(labels.lbl_mer_gwa).mean()
        wp_loc_mer_gwa = xr.open_mfdataset(results_pathg+"/windpower_??_MERRA2_GWA.nc")
        shares = xr.DataArray(BPA_mer_gwa.p[BPA_mer_gwa.p>0].values,dims = 'location',coords = {'location':range(sum(BPA_mer_gwa.p>0))})
        wp_BPA_mer_gwa = (wp_loc_mer_gwa.isel(location = np.where(BPA_mer_gwa.p.values > 0)[0]).assign_coords(location = range(len(shares))) * shares).sum('location')
        wp_BPA_mer_gwa.to_netcdf(ofile)
        del(BPA_mer_gwa,wp_loc_mer_gwa,wp_BPA_mer_gwa)
    
    # ERA5
    ofile = results_path + "/windpower_BPA_ERA5.nc"
    if ofile not in glob.glob(results_path + "/*"):
        BPA_era = pBPA.groupby(labels.lbl_era).mean()
        wp_loc_era = xr.open_dataset(results_path+"/windpower_stat_ERA5.nc")
        wp_BPA_era = (wp_loc_era * np.array(BPA_era.p)).sum('location')
        wp_BPA_era.to_netcdf(ofile)
        del(BPA_era,wp_loc_era,wp_BPA_era)
    
    # ERA5 + GWA
    ofile = results_pathg + "/windpower_BPA_ERA5_GWA.nc"
    if ofile not in glob.glob(results_pathg + "/*"):
        BPA_era_gwa = pBPA.groupby(labels.lbl_era_gwa).mean()
        wp_loc_era_gwa = xr.open_mfdataset(results_pathg+"/windpower_??_ERA5_GWA.nc")
        shares = xr.DataArray(BPA_era_gwa.p[BPA_era_gwa.p>0].values,dims = 'location',coords = {'location':range(sum(BPA_era_gwa.p>0))})
        wp_BPA_era_gwa = (wp_loc_era_gwa.isel(location = np.where(BPA_era_gwa.p.values > 0)[0]).assign_coords(location = range(len(shares))) * shares).sum('location')
        wp_BPA_era_gwa.to_netcdf(ofile)
        del(BPA_era_gwa,wp_loc_era_gwa,wp_BPA_era_gwa)
        
    del(BPA_parks,windturbines,labels,pBPA)


# sum up New England
if len(glob.glob(results_pathg + "/windpower_NewEngland_*.nc")) < 4:
    print("summing up New England...")
    NE_states = ['CT','NH','ME','MA','RI','VT']

    # MERRA2
    ofile = results_path+"/windpower_NewEngland_MERRA2.nc"
    if ofile not in glob.glob(results_path + "/*"):
        wp_state = xr.open_dataset(results_path+"/windpower_states_MERRA2.nc")
        wpNE = wp_state.sel(state=NE_states).sum('state')
        wpNE.to_netcdf(ofile)
        del(wp_state,wpNE)

    # MERRA2 + GWA
    ofile = results_pathg+"/windpower_NewEngland_MERRA2_GWA.nc"
    if ofile not in glob.glob(results_pathg + "/*"):
        wp_state = xr.open_dataset(results_pathg+"/windpower_states_MERRA2_GWA.nc")
        wpNE = wp_state.sel(state=NE_states).sum('state')
        wpNE.to_netcdf(ofile)
        del(wp_state,wpNE)

    # ERA5
    ofile = results_path+"/windpower_NewEngland_ERA5.nc"
    if ofile not in glob.glob(results_path + "/*"):
        wp_state = xr.open_dataset(results_path+"/windpower_states_ERA5.nc")
        wpNE = wp_state.sel(state=NE_states).sum('state')
        wpNE.to_netcdf(ofile)
        del(wp_state,wpNE)

    # ERA5 + GWA
    ofile = results_pathg+"/windpower_NewEngland_ERA5_GWA.nc"
    if ofile not in glob.glob(results_pathg + "/*"):
        wp_state = xr.open_dataset(results_pathg+"/windpower_states_ERA5_GWA.nc")
        wpNE = wp_state.sel(state=NE_states).sum('state')
        wpNE.to_netcdf(ofile)
        del(wp_state,wpNE)


# sum up USA
if len(glob.glob(results_pathg + "/windpower_USA_*.nc")) < 4:
    print("summing up USA...")

    # MERRA2
    ofile = results_path+"/windpower_USA_MERRA2.nc"
    if ofile not in glob.glob(results_path + "/*"):
        wp_state = xr.open_dataset(results_path+"/windpower_states_MERRA2.nc")
        wpUSA = wp_state.sum('state')
        wpUSA.to_netcdf(ofile)

    # MERRA2 + GWA
    ofile = results_pathg+"/windpower_USA_MERRA2_GWA.nc"
    if ofile not in glob.glob(results_pathg + "/*"):
        wp_state = xr.open_dataset(results_pathg+"/windpower_states_MERRA2_GWA.nc")
        wpUSA = wp_state.sum('state')
        wpUSA.to_netcdf(ofile)

    # ERA5
    ofile = results_path+"/windpower_USA_ERA5.nc"
    if ofile not in glob.glob(results_path + "/*"):
        wp_state = xr.open_dataset(results_path+"/windpower_states_ERA5.nc")
        wpUSA = wp_state.sum('state')
        wpUSA.to_netcdf(ofile)

    # ERA5 + GWA
    ofile = results_pathg+"/windpower_USA_ERA5_GWA.nc"
    if ofile not in glob.glob(results_pathg + "/*"):
        wp_state = xr.open_dataset(results_pathg+"/windpower_states_ERA5_GWA.nc")
        wpUSA = wp_state.sum('state')
        wpUSA.to_netcdf(ofile)

