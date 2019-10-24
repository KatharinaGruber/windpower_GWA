import math
import numpy as np
import pandas as pd
import statsmodels.api as sm
import xarray as xr

from functools import reduce

from paths_usa import *

windturbines = pd.read_csv(usa_path+"/windturbines_usa.csv",delimiter=';')
# fill in missing data
# commissioning dates
y = windturbines.p_year.fillna(windturbines.p_year.mean().round(0))
# capacities - fill with yearly mean
wt = windturbines.dropna(subset=['t_cap'])
y_cap = wt.t_cap.groupby(wt.p_year).mean().reindex(range(1982,2020)).ffill()
cap = windturbines.set_index(y.values).t_cap.fillna(y_cap).values
cap[np.isnan(cap)] = y_cap[1982] # fill missing values of 1981 with value of 1982, because none for 1981
# hub heights - fill with yearly mean
wt = windturbines.dropna(subset=['t_hh'])
y_hh = wt.t_hh.groupby(wt.p_year).mean().reindex(range(1981,2020)).ffill().bfill()
hh = windturbines.set_index(y.values).t_hh.fillna(y_hh).values

# make commissioning dates - year is given -> use middle of year
t = [np.datetime64(str(int(year))+"-06-01T00:00:00") for year in y]

# calculate specific power - regress rotor diameter from height where both are given
hh_rd = windturbines[['t_hh','t_rd']].dropna()
model = sm.OLS(hh_rd.t_rd, hh_rd.t_hh).fit()
mis_rd = model.predict(hh[np.isnan(windturbines.t_rd)])
rd = windturbines.t_rd.copy(deep=True)
rd[np.isnan(windturbines.t_rd)] = mis_rd
sp = cap/(rd**2*math.pi/4)*1000

# determine the capacities of turbines with sp < 100 (powercurve does not look realistic with this value)
caps_newsp = np.sort(np.unique(cap[sp<100]))
# some are not in rest of data - find closest value
diff = abs(cap[sp>=100]-caps_newsp[[cn not in cap[sp>=100] for cn in caps_newsp]])
caps_newsp[[cn not in cap[sp>=100] for cn in caps_newsp]] = cap[sp>=100][np.where(diff==min(diff))][0]
# get new sp
new_sp = [sp[(sp>=100)&(cap==c)].mean() for c in caps_newsp]
caps_newsp = np.sort(np.unique(cap[sp<100]))
# insert new sp values
sp[sp<100] = pd.Series(cap)[sp<100].map(dict(zip(caps_newsp,new_sp)))

# join data making a multiindex with time and state
time_state = tuples = tuple(zip(list(windturbines.t_state),t))
mi = pd.MultiIndex.from_tuples(time_state, names=['state', 'time'])
turbine_data = pd.DataFrame({'capacity':cap,
                             'height':hh,
                             'lon':windturbines.xlong.values,
                             'lat':windturbines.ylat.values,
                             'sp':sp.values},index = mi)

# remove Guam because it's outside
turbine_data = turbine_data.drop('GU',axis=0,level=0)

# add index
turbine_data['ind'] = range(turbine_data.shape[0])


# MERRA data preparation
##########################
data_mer = xr.open_dataset(mer_path+"/merra2_wind_USA_200012.nc")
# Create dataframe with sequence the size of MERRA-2 grid to find out which turbines interpolate to the same point
in_seq = xr.Dataset({'x':(['lat','lon'],
                          np.array(range(data_mer.DISPH.isel(time=0).values.size)).reshape(data_mer.DISPH.isel(time=0).values.shape))},
                     coords = {'lat':data_mer.lat.values,
                               'lon':data_mer.lon.values})
# interpolate to indices
ip = in_seq.interp(coords={"lon":xr.DataArray(turbine_data.lon,dims='location'),
                           "lat":xr.DataArray(turbine_data.lat,dims='location')},method="nearest").to_dataframe()

# Load GWA data and extract values at locations to find unique locations					   
GWA_USA = xr.open_rasterio(usa_path+'/GWA/GWA_USA50m.tif')
# three regions are not included in the continental USA region: Alaska, Haitii and Puerto Rico
GWA_AK = xr.open_rasterio(usa_path+'/GWA/GWA_AK50m.tif')
GWA_HI = xr.open_rasterio(usa_path+'/GWA/GWA_HI50m.tif')
GWA_PR = xr.open_rasterio(usa_path+'/GWA/GWA_PR50m.tif')
# interpolate GWA to all locations and compare if same locations from rea dataset match same locations from GWA
GWA_ip = GWA_USA.sel(band=1).interp(coords={"x":xr.DataArray(turbine_data.lon,dims='location'),
                                            "y":xr.DataArray(turbine_data.lat,dims='location')},
                                    method="nearest").to_dataframe(name='GWA')
GWA_ak = GWA_AK.sel(band=1).interp(coords={"x":xr.DataArray(turbine_data.xs('AK').lon,dims='location'),
                                           "y":xr.DataArray(turbine_data.xs('AK').lat,dims='location')},
                                   method="nearest").to_dataframe(name='GWA')
GWA_ip.loc[GWA_ip.index.get_level_values(0).values=='AK','GWA'] = GWA_ak.GWA.values
GWA_hi = GWA_HI.sel(band=1).interp(coords={"x":xr.DataArray(turbine_data.xs('HI').lon,dims='location'),
                                           "y":xr.DataArray(turbine_data.xs('HI').lat,dims='location')},
                                   method="nearest").to_dataframe(name='GWA')
GWA_ip.loc[GWA_ip.index.get_level_values(0).values=='HI','GWA'] = GWA_hi.GWA.values
GWA_pr = GWA_PR.sel(band=1).interp(coords={"x":xr.DataArray(turbine_data.xs('PR').lon,dims='location'),
                                           "y":xr.DataArray(turbine_data.xs('PR').lat,dims='location')},
                                   method="nearest").to_dataframe(name='GWA')
GWA_ip.loc[GWA_ip.index.get_level_values(0).values=='PR','GWA'] = GWA_pr.GWA.values
# create labels for accumulating data to determine which configuartions are unique label consists of:
# - state (because some windturbines interpolate to same point in grid but are in different states)
# - interpolation point
# - commissioning date
# - specific power
# - hub height
lbl_mer = reduce(np.core.defchararray.add,[np.array(ip.index.get_level_values(0)).astype('U'), "_",
                                           np.array(ip.x).astype('U'),'_',
                                           ip.index.get_level_values(1).astype('U'),"_sp",
                                           turbine_data.sp.astype('U'),"_hh",
                                           turbine_data.height.astype('U')])
lbl_mer_gwa = reduce(np.core.defchararray.add,[np.array(ip.index.get_level_values(0)).astype('U'), "_",
                                               np.array(ip.x).astype('U'),'_',
                                               ip.index.get_level_values(1).astype('U'),"_sp",
                                               turbine_data.sp.astype('U'),"_hh",
                                               turbine_data.height.astype('U'), "_gwa",
                                               GWA_ip.GWA.astype('U')])

# Aggregate turbine data - sum up capacities for same locations and same specifications (hubheight, specific power)
# without GWA
cap_ag = turbine_data.capacity.groupby(lbl_mer).sum().values
ind_ag = turbine_data.ind.groupby(lbl_mer).min().values
state_ag = turbine_data.index.get_level_values(0).values[ind_ag]
comm_ag = turbine_data.index.get_level_values(1).values[ind_ag]
lon_ag = turbine_data.lon.values[ind_ag]
lat_ag = turbine_data.lat.values[ind_ag]
sp_ag = turbine_data.sp.values[ind_ag]
hh_ag = turbine_data.height.values[ind_ag]
turbine_data_mer = pd.DataFrame({'state':state_ag,
                                 'commissioning':comm_ag,
                                 'capacity':cap_ag,
                                 'lon':lon_ag,
                                 'lat':lat_ag,
                                 'sp':sp_ag,
                                 'height':hh_ag})
# with GWA
cap_ag = turbine_data.capacity.groupby(lbl_mer_gwa).sum().values
ind_ag = turbine_data.ind.groupby(lbl_mer_gwa).min().values
state_ag = turbine_data.index.get_level_values(0).values[ind_ag]
comm_ag = turbine_data.index.get_level_values(1).values[ind_ag]
lon_ag = turbine_data.lon.values[ind_ag]
lat_ag = turbine_data.lat.values[ind_ag]
sp_ag = turbine_data.sp.values[ind_ag]
hh_ag = turbine_data.height.values[ind_ag]
turbine_data_mer_gwa = pd.DataFrame({'state':state_ag,
                                     'commissioning':comm_ag,
                                     'capacity':cap_ag,
                                     'lon':lon_ag,
                                     'lat':lat_ag,
                                     'sp':sp_ag,
                                     'height':hh_ag})


# ERA data preparation
##########################
data_era = xr.open_dataset(era_path+"/era5_wind_USA_200012.nc")
# Create dataframe with sequence the size of ERA5 grid to find out which turbines interpolate to the same point
in_seq = xr.Dataset({'x':(['lat','lon'],
                          np.array(range(data_era.u10.isel(time=0).values.size)).reshape(data_era.u10.isel(time=0).values.shape))},
                     coords = {'lat':data_era.latitude.values,
                               'lon':data_era.longitude.values})
# interpolate to indices
ip = in_seq.interp(coords={"lon":xr.DataArray(turbine_data.lon,dims='location'),
                           "lat":xr.DataArray(turbine_data.lat,dims='location')},method="nearest").to_dataframe()

# Load GWA data and extract values at locations to find unique locations					   
GWA_USA = xr.open_rasterio(usa_path+'/GWA/GWA_USA100m.tif')
# three regions are not included in the continental USA region: Alaska, Haitii and Puerto Rico
GWA_AK = xr.open_rasterio(usa_path+'/GWA/GWA_AK100m.tif')
GWA_HI = xr.open_rasterio(usa_path+'/GWA/GWA_HI100m.tif')
GWA_PR = xr.open_rasterio(usa_path+'/GWA/GWA_PR100m.tif')
# interpolate GWA to all locations and compare if same locations from rea dataset match same locations from GWA
GWA_ip = GWA_USA.sel(band=1).interp(coords={"x":xr.DataArray(turbine_data.lon,dims='location'),
                                            "y":xr.DataArray(turbine_data.lat,dims='location')},
                                    method="nearest").to_dataframe(name='GWA')
GWA_ak = GWA_AK.sel(band=1).interp(coords={"x":xr.DataArray(turbine_data.xs('AK').lon,dims='location'),
                                           "y":xr.DataArray(turbine_data.xs('AK').lat,dims='location')},
                                   method="nearest").to_dataframe(name='GWA')
GWA_ip.loc[GWA_ip.index.get_level_values(0).values=='AK','GWA'] = GWA_ak.GWA.values
GWA_hi = GWA_HI.sel(band=1).interp(coords={"x":xr.DataArray(turbine_data.xs('HI').lon,dims='location'),
                                           "y":xr.DataArray(turbine_data.xs('HI').lat,dims='location')},
                                   method="nearest").to_dataframe(name='GWA')
GWA_ip.loc[GWA_ip.index.get_level_values(0).values=='HI','GWA'] = GWA_hi.GWA.values
GWA_pr = GWA_PR.sel(band=1).interp(coords={"x":xr.DataArray(turbine_data.xs('PR').lon,dims='location'),
                                           "y":xr.DataArray(turbine_data.xs('PR').lat,dims='location')},
                                   method="nearest").to_dataframe(name='GWA')
GWA_ip.loc[GWA_ip.index.get_level_values(0).values=='PR','GWA'] = GWA_pr.GWA.values
# create labels for accumulating data to determine which configuartions are unique label consists of:
# - state (because some windturbines interpolate to same point in grid but are in different states)
# - interpolation point
# - commissioning date
# - specific power
# - hub height
lbl_era = reduce(np.core.defchararray.add,[np.array(ip.index.get_level_values(0)).astype('U'), "_",
                                           np.array(ip.x).astype('U'),'_',
                                           ip.index.get_level_values(1).astype('U'),"_sp",
                                           turbine_data.sp.astype('U'),"_hh",
                                           turbine_data.height.astype('U')])
lbl_era_gwa = reduce(np.core.defchararray.add,[np.array(ip.index.get_level_values(0)).astype('U'), "_",
                                               np.array(ip.x).astype('U'),'_',
                                               ip.index.get_level_values(1).astype('U'),"_sp",
                                               turbine_data.sp.astype('U'),"_hh",
                                               turbine_data.height.astype('U'), "_gwa",
                                               GWA_ip.GWA.astype('U')])

# Aggregate turbine data - sum up capacities for same locations and same specifications (hubheight, specific power)
# without GWA
cap_ag = turbine_data.capacity.groupby(lbl_era).sum().values
ind_ag = turbine_data.ind.groupby(lbl_era).min().values
state_ag = turbine_data.index.get_level_values(0).values[ind_ag]
comm_ag = turbine_data.index.get_level_values(1).values[ind_ag]
lon_ag = turbine_data.lon.values[ind_ag]
lat_ag = turbine_data.lat.values[ind_ag]
sp_ag = turbine_data.sp.values[ind_ag]
hh_ag = turbine_data.height.values[ind_ag]
turbine_data_era = pd.DataFrame({'state':state_ag,
                                 'commissioning':comm_ag,
                                 'capacity':cap_ag,
                                 'lon':lon_ag,
                                 'lat':lat_ag,
                                 'sp':sp_ag,
                                 'height':hh_ag})
# with GWA
cap_ag = turbine_data.capacity.groupby(lbl_era_gwa).sum().values
ind_ag = turbine_data.ind.groupby(lbl_era_gwa).min().values
state_ag = turbine_data.index.get_level_values(0).values[ind_ag]
comm_ag = turbine_data.index.get_level_values(1).values[ind_ag]
lon_ag = turbine_data.lon.values[ind_ag]
lat_ag = turbine_data.lat.values[ind_ag]
sp_ag = turbine_data.sp.values[ind_ag]
hh_ag = turbine_data.height.values[ind_ag]
turbine_data_era_gwa = pd.DataFrame({'state':state_ag,
                                     'commissioning':comm_ag,
                                     'capacity':cap_ag,
                                     'lon':lon_ag,
                                     'lat':lat_ag,
                                     'sp':sp_ag,
                                     'height':hh_ag})
									 
# save to files
turbine_data_mer.to_csv(usa_path + '/turbine_data_mer.csv')
turbine_data_mer_gwa.to_csv(usa_path + '/turbine_data_mer_gwa.csv')
turbine_data_era.to_csv(usa_path + '/turbine_data_era.csv')
turbine_data_era_gwa.to_csv(usa_path + '/turbine_data_era_gwa.csv')

# merge labels and save
labels = pd.DataFrame({'lbl_mer':lbl_mer,
                       'lbl_mer_gwa':lbl_mer_gwa,
                       'lbl_era':lbl_era,
                       'lbl_era_gwa':lbl_era_gwa})
labels.to_csv(usa_path + '/labels_turbine_data.csv')