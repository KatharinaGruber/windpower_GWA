#!/usr/bin/env python
# coding: utf-8
"""Wind power simulation in South Africa with Global Wind Atlas bias correction

This script simulates wind power generation in South Africa
from ERA5 reanalysis data and performs wind speed bias correction
with DTU's Global Wind Atlas.

First, ERA5 data are downloaded, for this the CDS-API must be installed:
https://cds.climate.copernicus.eu/api-how-to

For simulating wind power generation, wind park data is needed.
The data have been aggregated by several sources, including:
REDIS: http://redis.energy.gov.za/power-producers/
Global Power Plant Database:
    http://datasets.wri.org/dataset/globalpowerplantdatabase
Wikipedia: https://en.wikipedia.org/wiki/List_of_wind_farms_in_South_Africa

Download Global Wind Atlas wind speeds at 100 m height for South Africa here:
    https://globalwindatlas.info/

Downlaod wind power generation data from REDIS:
    http://redis.energy.gov.za/power-production/

"""

import os

import numpy as np
import pandas as pd
import xarray as xr

from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

# define path for windpark data
windparks_path = "C:/..."
# define path for Global Wind Atlas (GWA) data
GWA_path = "C:/..."
# define ERA5 path
era5_path = "C:/..."
# define path of wind power generation data
windpower_path = "C:/..."
# define path for results
results_path = "C:/..."


# download ERA5 data
os.chdir(era5_path)
os.system("era5_download_zaf.py")


# define chunks for loading ERA5 data
# only chunk over time, because cannot chunk over dimension where interpolated
chk = {
    "time": 36
}
# load ERA5 data
data = xr.open_mfdataset(era5_path+"/*.nc", chunks=chk)


# load locations of wind power plants
windparks = pd.read_csv(windparks_path+"/windparks_southafrica.csv",
                        delimiter=';')

# extract years, months and days of commissioning dates
# if month not given set to 6, if day not given set to 15
y = np.array([x[:4] for x in windparks.commissioning])
m = np.array([x[5:7] for x in windparks.commissioning])
m[m == ''] = '06'
d = np.array([x[8:10] for x in windparks.commissioning])
d[d == ''] = '15'
# join years, months and days to datetime format
t = [np.datetime64(y[i] + '-' + m[i] + '-' + d[i] + "T00:00:00")
     for i in np.arange(len(y))]

# extract longitudes and latitudes of locations of wind parks
lons = windparks.Longitude
lats = windparks.Latitude


# create multiindex of locations of windparks with 2 indices:
# longitude and latitude
idx = pd.MultiIndex.from_arrays([lons,
                                 lats],
                                names=['lon', 'lat'])


# interpolate wind speeds to locations of wind parks
# store in list and concatenate to xarray dataset
# (for some reason using the same method as for the GWA
# results in a warning and may produce errors in the future)
list_wind_wp = []
for i in np.arange(len(lons)):
    wind_1wp = data.interp(coords={"longitude": lons[i],
                                   "latitude": lats[i]}, method="nearest")
    list_wind_wp.append(wind_1wp)
wind_windparks = xr.concat(list_wind_wp, dim=idx)

# remove unnecessary duplicate longitude and latitude indices
wind_windparks = wind_windparks.drop('longitude').drop('latitude')
# and rename concat_dim to location
wind_windparks = wind_windparks.rename({'concat_dim': 'location'})


# -------------------------------------------------------------
# Inter- and extrapolation of wind speeds to hub height (108 m)
# -------------------------------------------------------------

# calculate wind at 10 and 100m height with pythagoras
# from u and v wind speed component:
# v_{eff} = \sqrt{u^2+v^2}

# wind in 10m height calculated from u and v component
windh10 = np.sqrt(wind_windparks.v10.values**2
                  + wind_windparks.u10.values**2)
# wind in 100m height calculated from u and v component
windh100 = np.sqrt(wind_windparks.v100.values**2
                   + wind_windparks.u100.values**2)


# calculate alpha friction coefficient
# \alpha = \frac{(ln(v_2/v_1)}{(ln(h_2/h_1))}
alpha = (np.log(windh100/windh10))/(np.log(100/10))


# then calculate wind speed at height 108m (height of Enercon E-82)
# v_2 = v_1 * \Bigl(\frac{h_2}{h_1}\Bigr)^\alpha
turbine_height = 108
windh108 = windh10 * (turbine_height/10)**alpha

# -------------------------------------------------------------
# calculate wind power generation from wind speeds
# -------------------------------------------------------------
# Enercon E-82 power curve
# https://www.enercon.de/fileadmin/Redakteur/Medien-Portal/broschueren/pdf/en/ENERCON_Produkt_en_06_2015.pdf
wind_speeds = (np.arange(0, 26, step=1.0))
generation_kw = [0.0, 0.000000000001, 3.0, 25.0, 82.0, 175.0, 321.0, 532.0,
                 815.0, 1180.0, 1580.0, 1810.0, 1980.0] + 13 * [2050.0]
plt.plot(wind_speeds, generation_kw)

# create power curve function,
# which linearly interpolates between power curve points
power_curve = interp1d(wind_speeds, generation_kw)

# calculate wind power generation
wp1 = xr.apply_ufunc(power_curve, windh108,
                     dask='parallelized',
                     output_dtypes=[np.float64])
# fetch installed capacity and divide by 2000
# to make factor for capacity of Enercon E-82
cap = list(windparks.Capacity/2000.0)
# multiply with installed capacity
wp2 = cap*wp1.transpose()


def commission(date1, timeseries, windpower):
    """ commission wind parks
    function which "commissions" wind park :
    before commissioning production is set to 0
    parameters:
    - date1: commissioning date, dtype datetime64
    - timeseries: time series of period, dtype datetime64
    - windpower: wind power scaled by installed capacity,
                 length same as timeseries
    """
    bin_sel = np.zeros(len(timeseries))
    bin_sel[timeseries >= date1] = 1
    out = windpower * bin_sel
    return out


# prepare array for commissioned wind power
wp3 = np.zeros(wp2.shape)
# make wind power generation start at commissioning date
for i in np.arange(len(t)):
    wp3[:, i] = commission(t[i], wind_windparks.time.values, wp2[:, i])


# make into xarray dataset
windpower_windparks = xr.DataArray(wp3,
                                   coords={'time': wind_windparks.time.values,
                                           'location': idx},
                                   dims=('time',
                                         'location')
                                   ).to_dataset(name='wind_power')


# -------------------------------------------------------------
# Preparation for statistical analysis
# -------------------------------------------------------------

# extract the indices of locations from different areas
ind_nc = np.where(windparks.Area == 'Northern Cape')
ind_ec = np.where(windparks.Area == 'Eastern Cape')
ind_wc = np.where(windparks.Area == 'Western Cape')


# sum the wind power generation for these regions
windpower_nc = windpower_windparks.isel(location=ind_nc[0].tolist()).wind_power.values
windpower_ec = windpower_windparks.isel(location=ind_ec[0].tolist()).wind_power.values
windpower_wc = windpower_windparks.isel(location=ind_wc[0].tolist()).wind_power.values
windpower_kW_area = pd.DataFrame({'Northern_Cape': np.sum(windpower_nc, axis=1),
                                  'Eastern_Cape': np.sum(windpower_ec, axis=1),
                                  'Western_Cape': np.sum(windpower_wc, axis=1)},
                                 index=windpower_windparks.time.values)

# Load data Eastern Cape
data_ec = []
for year in range(2014, 2019):
    datay = pd.read_csv(windpower_path
                        + "/EasternCape/Hourly_Electricity_production_[Load_Factor_[[%]]_data"
                        + str(year) + ".csv", delimiter=';')
    if(len(data_ec)):
        data_ec = pd.concat([data_ec, datay])
    else:
        data_ec = datay

# Load data Western Cape
data_wc = []
for year in range(2015, 2019):
    datay = pd.read_csv(windpower_path
                        + "/WesternCape/Hourly_Electricity_production_[Load_Factor_[[%]]_data"
                        + str(year) + ".csv", delimiter=';')
    if(len(data_wc)):
        data_wc = pd.concat([data_wc, datay])
    else:
        data_wc = datay

# Load data Northern Cape
data_nc = []
for year in range(2017, 2019):
    datay = pd.read_csv(windpower_path
                        + "/NorthernCape/Hourly_Electricity_production_[Load_Factor_[[%]]_data"
                        + str(year) + ".csv", delimiter=';')
    if(len(data_nc)):
        data_nc = pd.concat([data_nc, datay])
    else:
        data_nc = datay


# -------------------------------------------------------------
# Arrange resulting time series
# -------------------------------------------------------------

# create pandas dataframes for each of the locations

# Eastern Cape
y = np.array([data_ec['Settlement DateTime'].iloc[x][6:10]
              for x in range(data_ec.shape[0])])
m = np.array([data_ec['Settlement DateTime'].iloc[x][3:5]
              for x in range(data_ec.shape[0])])
d = np.array([data_ec['Settlement DateTime'].iloc[x][:2]
              for x in range(data_ec.shape[0])])
h = np.array([data_ec['Settlement DateTime'].iloc[x][11:13]
              for x in range(data_ec.shape[0])])
t1 = [np.datetime64(y[i]+'-'+m[i]+'-'+d[i]+"T"+h[i]+":00:00")
      for i in np.arange(len(y))]
# create dataframe with datetime index and simulated data
# and create column to fill in historical data
production_MW_ec = pd.DataFrame({'wp_MWh': data_ec['Production MWh'].tolist(),
                                 'sim_wp_MWh': np.nan},
                                index=t1)

# Western Cape
y = np.array([data_wc['Settlement DateTime'].iloc[x][6:10]
              for x in range(data_wc.shape[0])])
m = np.array([data_wc['Settlement DateTime'].iloc[x][3:5]
              for x in range(data_wc.shape[0])])
d = np.array([data_wc['Settlement DateTime'].iloc[x][:2]
              for x in range(data_wc.shape[0])])
h = np.array([data_wc['Settlement DateTime'].iloc[x][11:13]
              for x in range(data_wc.shape[0])])
t1 = [np.datetime64(y[i]+'-'+m[i]+'-'+d[i]+"T"+h[i]+":00:00")
      for i in np.arange(len(y))]
# create dataframe with datetime index and simulated data
# and create column to fill in historical data
production_MW_wc = pd.DataFrame({'wp_MWh': data_wc['Production MWh'].tolist(),
                                 'sim_wp_MWh': np.nan},
                                index=t1)

# Northern Cape
y = np.array([data_nc['Settlement DateTime'].iloc[x][6:10]
              for x in range(data_nc.shape[0])])
m = np.array([data_nc['Settlement DateTime'].iloc[x][3:5]
              for x in range(data_nc.shape[0])])
d = np.array([data_nc['Settlement DateTime'].iloc[x][:2]
              for x in range(data_nc.shape[0])])
h = np.array([data_nc['Settlement DateTime'].iloc[x][11:13]
              for x in range(data_nc.shape[0])])
t1 = [np.datetime64(y[i]+'-'+m[i]+'-'+d[i]+"T"+h[i]+":00:00")
      for i in np.arange(len(y))]
# create dataframe with datetime index and simulated data
# and create column to fill in historical data
production_MW_nc = pd.DataFrame({'wp_MWh': data_nc['Production MWh'].tolist(),
                                 'sim_wp_MWh': np.nan},
                                index=t1)

# fill in simulated data - cannot be added simply as other column,
# as in historical data some time steps are missing
# therefore match simulated data to time steps
production_MW_ec.sim_wp_MWh = production_MW_ec.index.map(windpower_kW_area.Eastern_Cape)
production_MW_wc.sim_wp_MWh = production_MW_wc.index.map(windpower_kW_area.Western_Cape)
production_MW_nc.sim_wp_MWh = production_MW_nc.index.map(windpower_kW_area.Northern_Cape)

# remove missing values
production_MW_ec = production_MW_ec.dropna(axis=0)
production_MW_wc = production_MW_wc.dropna(axis=0)
production_MW_nc = production_MW_nc.dropna(axis=0)


# -------------------------------------------------------------
# Global Wind Atlas bias correction
# -------------------------------------------------------------

# Get Global Wind Atlas data
# use rasterio, it's the easiest way and directly
# makes it into an xarray for easy interpolation
GWA2 = xr.open_rasterio(GWA_path+'/GWA_ZAF.tif')

# interpolate to locations of wind turbines
GWA_locations = GWA2.interp(coords={"x": xr.DataArray(lons, dims='location'),
                                    "y": xr.DataArray(lats, dims='location')},
                            method="linear")

# calculate correction factor wind speeds by division
corr_fac_GWA = GWA_locations.values/windh100.mean(axis=1)

# apply correction factor
windh108_GWA = windh108.transpose() * corr_fac_GWA

# calculate wind power generation like before
windh108_GWAc = windh108_GWA
windh108_GWAc[np.where(windh108_GWA > 25)] = 25
wp1_GWA = xr.apply_ufunc(power_curve,
                         windh108_GWAc,
                         dask='parallelized',
                         output_dtypes=[np.float64])

# multiply with installed capacity
wp2_GWA = cap*wp1_GWA

# make wind power generation start at commissioning date
wp3_GWA = np.zeros(wp2_GWA.shape)
for i in np.arange(len(t)):
    wp3_GWA[:, i] = commission(t[i], wind_windparks.time.values, wp2_GWA[:, i])

# create xarray dataset
windpower_windparks_GWA = xr.DataArray(wp3_GWA,
                                       coords={'time': wind_windparks.time.values,
                                               'location': idx},
                                       dims=('time',
                                             'location')
                                       ).to_dataset(name='wind_power')

# -------------------------------------------------------------
# Preparation of GWA simulation results for statistical analysis
# -------------------------------------------------------------

# sum the wind power generation for these regions
windpower_GWA_nc = windpower_windparks_GWA.isel(location=ind_nc[0].tolist()).wind_power.values
windpower_GWA_ec = windpower_windparks_GWA.isel(location=ind_ec[0].tolist()).wind_power.values
windpower_GWA_wc = windpower_windparks_GWA.isel(location=ind_wc[0].tolist()).wind_power.values
windpower_kW_area_GWA = pd.DataFrame({'Northern_Cape': np.sum(windpower_GWA_nc, axis=1),
                                      'Eastern_Cape': np.sum(windpower_GWA_nc, axis=1),
                                      'Western_Cape': np.sum(windpower_GWA_nc, axis=1)},
                                     index=windpower_windparks_GWA.time.values)


# copy prepared dataframe from previous analyses
# deepcopy because otherwise it will overwrite previous results
production_MW_ec_GWA = production_MW_ec.copy(deep=True)
production_MW_wc_GWA = production_MW_wc.copy(deep=True)
production_MW_nc_GWA = production_MW_nc.copy(deep=True)
# fill in simulated data - cannot be added simply as other column,
# as in historical data some time steps are missing
# therefore match simulated data to time steps
production_MW_ec_GWA.sim_wp_MWh = production_MW_ec_GWA.index.map(windpower_kW_area_GWA.Eastern_Cape)
production_MW_wc_GWA.sim_wp_MWh = production_MW_wc_GWA.index.map(windpower_kW_area_GWA.Western_Cape)
production_MW_nc_GWA.sim_wp_MWh = production_MW_nc_GWA.index.map(windpower_kW_area_GWA.Northern_Cape)
# remove missing values
production_MW_ec_GWA = production_MW_ec_GWA.dropna(axis=0)
production_MW_wc_GWA = production_MW_wc_GWA.dropna(axis=0)
production_MW_nc_GWA = production_MW_nc_GWA.dropna(axis=0)


# -------------------------------------------------------------
# Statistical analysis for regions (capes)
# -------------------------------------------------------------
cors = [np.corrcoef(production_MW_ec.wp_MWh,
                    production_MW_ec.sim_wp_MWh)[0, 1],
        np.corrcoef(production_MW_wc.wp_MWh,
                    production_MW_wc.sim_wp_MWh)[0, 1],
        np.corrcoef(production_MW_nc.wp_MWh,
                    production_MW_nc.sim_wp_MWh)[0, 1]]

rmses = [np.sqrt(np.mean((np.array(production_MW_ec.wp_MWh)
                          - np.array(production_MW_ec.sim_wp_MWh))**2)),
         np.sqrt(np.mean((np.array(production_MW_wc.wp_MWh)
                          - np.array(production_MW_wc.sim_wp_MWh))**2)),
         np.sqrt(np.mean((np.array(production_MW_nc.wp_MWh)
                          - np.array(production_MW_nc.sim_wp_MWh))**2))]

mbes = [(np.array(production_MW_ec.wp_MWh)
         - np.array(production_MW_ec.sim_wp_MWh)).mean(),
        (np.array(production_MW_wc.wp_MWh)
         - np.array(production_MW_wc.sim_wp_MWh)).mean(),
        (np.array(production_MW_nc.wp_MWh)
        - np.array(production_MW_nc.sim_wp_MWh)).mean()]

means = [[production_MW_ec.wp_MWh.mean(), production_MW_ec.sim_wp_MWh.mean()],
         [production_MW_wc.wp_MWh.mean(), production_MW_wc.sim_wp_MWh.mean()],
         [production_MW_nc.wp_MWh.mean(), production_MW_nc.sim_wp_MWh.mean()]]


# with Global Wind Atlas
corsGWA = [np.corrcoef(production_MW_ec_GWA.wp_MWh,
                       production_MW_ec_GWA.sim_wp_MWh)[0, 1],
           np.corrcoef(production_MW_wc_GWA.wp_MWh,
                       production_MW_wc_GWA.sim_wp_MWh)[0, 1],
           np.corrcoef(production_MW_nc_GWA.wp_MWh,
                       production_MW_nc_GWA.sim_wp_MWh)[0, 1]]

rmsesGWA = [np.sqrt(np.mean((np.array(production_MW_ec_GWA.wp_MWh)
                             - np.array(production_MW_ec_GWA.sim_wp_MWh))**2)),
            np.sqrt(np.mean((np.array(production_MW_wc_GWA.wp_MWh)
                             - np.array(production_MW_wc_GWA.sim_wp_MWh))**2)),
            np.sqrt(np.mean((np.array(production_MW_nc_GWA.wp_MWh)
                             - np.array(production_MW_nc_GWA.sim_wp_MWh))**2))]

mbesGWA = [(np.array(production_MW_ec_GWA.wp_MWh)
            - np.array(production_MW_ec_GWA.sim_wp_MWh)).mean(),
           (np.array(production_MW_wc_GWA.wp_MWh)
            - np.array(production_MW_wc_GWA.sim_wp_MWh)).mean(),
           (np.array(production_MW_nc_GWA.wp_MWh)
            - np.array(production_MW_nc_GWA.sim_wp_MWh)).mean()]

meansGWA = [[production_MW_ec_GWA.wp_MWh.mean(),
             production_MW_ec_GWA.sim_wp_MWh.mean()],
            [production_MW_wc_GWA.wp_MWh.mean(),
             production_MW_wc_GWA.sim_wp_MWh.mean()],
            [production_MW_nc_GWA.wp_MWh.mean(),
             production_MW_nc_GWA.sim_wp_MWh.mean()]]


# merge the statistical results to dataframe and save as csv
pd.DataFrame({'Eastern_Cape': [cors[0], rmses[0], mbes[0], means[0][1]],
              'Eastern_Cape_GWA': [corsGWA[0], rmsesGWA[0], mbesGWA[0],
                                   meansGWA[0][1]],
              'Eastern_Cape_obs': [None] * 3 + [means[0][0]],
              'Western_Cape': [cors[1], rmses[1], mbes[1], means[1][1]],
              'Western_Cape_GWA': [corsGWA[1], rmsesGWA[1], mbesGWA[1],
                                   meansGWA[1][1]],
              'Western_Cape_obs': [None] * 3 + [means[1][0]],
              'Northern_Cape': [cors[2], rmses[2], mbes[2], means[2][1]],
              'Northern_Cape_GWA': [corsGWA[2], rmsesGWA[2], mbesGWA[2],
                                    meansGWA[2][1]],
              'Northern_Cape_obs': [None] * 3 + [means[2][0]]},
             index=['Correlation', 'RMSE_MWh', 'MBE_MWh', 'Mean_MWh'])

boxplot = [list(production_MW_ec.wp_MWh),
           list(production_MW_ec.sim_wp_MWh),
           list(production_MW_ec_GWA.sim_wp_MWh),
           list(production_MW_wc.wp_MWh),
           list(production_MW_wc.sim_wp_MWh),
           list(production_MW_wc_GWA.sim_wp_MWh),
           list(production_MW_nc.wp_MWh),
           list(production_MW_nc.sim_wp_MWh),
           list(production_MW_nc_GWA.sim_wp_MWh)]

fig = plt.figure(1, figsize=(12, 6))
# Create an axes instance
ax = fig.add_subplot(111)
# Create the boxplot
bp = ax.boxplot(boxplot)
ax.set_xticklabels(['EC-obs', 'EC-sim', 'EC-sim_GWA',
                    'WC-obs', 'WC-sim', 'WC-sim_GWA',
                    'NC-obs', 'NC-sim', 'NC-sim_GWA'])


# -------------------------------------------------------------
# Prepare data for statistical analysis for South Africa
# -------------------------------------------------------------

# need to remove days which are not complete, then aggregate daily
def drop_less_one_day(datetimeindex):
    '''
    function which creates series of 0 and 1
    when hours shall be kept (1, complete days with 24 hours)
    or hours shall be dropped (0, days where less than 24 hours are available)
    '''
    # datetimeindex to numbers
    dt = np.array(datetimeindex.year * 10000
                  + datetimeindex.month * 100
                  + datetimeindex.day)
    # count how often each day occurs
    length = np.unique(dt, return_counts=True)[1]
    # only keep days with 24 hours
    filt = np.array([0]*len(length))
    filt[np.equal(length, 24)] = 1
    # create code with 0 and 1 for hours
    # to keep (1) and discard (0)
    code = np.repeat(filt, length)

    return code


# ORDER before applying the drop less one day!
# for some reason there is a messup in the dates...
production_MW_ec = production_MW_ec.sort_index()
production_MW_wc = production_MW_wc.sort_index()
production_MW_nc = production_MW_nc.sort_index()
production_MW_ec_GWA = production_MW_ec_GWA.sort_index()
production_MW_wc_GWA = production_MW_wc_GWA.sort_index()
production_MW_nc_GWA = production_MW_nc_GWA.sort_index()

# examine which shares of data have been kept
print('Shares of data kept')
d = drop_less_one_day(production_MW_ec.index)
print("EC keep %: ", round(sum(np.equal(d, 1))
                           / len(production_MW_ec.index) * 100, 1))

d = drop_less_one_day(production_MW_wc.index)
print("WC keep %: ", round(sum(np.equal(d, 1))
                           / len(production_MW_wc.index) * 100, 1))

d = drop_less_one_day(production_MW_nc.index)
print("NC keep %: ", round(sum(np.equal(d, 1))
                           / len(production_MW_nc.index) * 100, 1))

# find out which lines to keep
production_MW_ec['keep'] = drop_less_one_day(production_MW_ec.index)
production_MW_ec_GWA['keep'] = drop_less_one_day(production_MW_ec_GWA.index)
production_MW_wc['keep'] = drop_less_one_day(production_MW_wc.index)
production_MW_wc_GWA['keep'] = drop_less_one_day(production_MW_wc_GWA.index)
production_MW_nc['keep'] = drop_less_one_day(production_MW_nc.index)
production_MW_nc_GWA['keep'] = drop_less_one_day(production_MW_nc_GWA.index)

# discard lines that do not belong to complete days
production_ec_filtered = production_MW_ec[production_MW_ec.keep == 1]
production_ecGWA_filtered = production_MW_ec_GWA[production_MW_ec_GWA.keep == 1]
production_wc_filtered = production_MW_wc[production_MW_wc.keep == 1]
production_wcGWA_filtered = production_MW_wc_GWA[production_MW_wc_GWA.keep == 1]
production_nc_filtered = production_MW_nc[production_MW_nc.keep == 1]
production_ncGWA_filtered = production_MW_nc_GWA[production_MW_nc_GWA.keep == 1]


# prepare dataframe for merging data from all regions and simulations
production_MW_all = pd.DataFrame({'EC_obs': np.nan,
                                  'EC_sim': windpower_kW_area.Eastern_Cape,
                                  'EC_simGWA': windpower_kW_area_GWA.Eastern_Cape,
                                  'WC_obs': np.nan,
                                  'WC_sim': windpower_kW_area.Western_Cape,
                                  'WC_simGWA': windpower_kW_area_GWA.Western_Cape,
                                  'NC_obs': np.nan,
                                  'NC_sim': windpower_kW_area.Northern_Cape,
                                  'NC_simGWA': windpower_kW_area_GWA.Northern_Cape, })


# fill in data values
production_MW_all.EC_obs = production_MW_all.index.map(production_ec_filtered.wp_MWh)
production_MW_all.WC_obs = production_MW_all.index.map(production_wc_filtered.wp_MWh)
production_MW_all.NC_obs = production_MW_all.index.map(production_nc_filtered.wp_MWh)


# fill in NANs at start with 0, because no production yet
production_MW_all.EC_obs[production_MW_all.index
                         < production_MW_all.EC_obs.first_valid_index()] = 0
production_MW_all.WC_obs[production_MW_all.index
                         < production_MW_all.WC_obs.first_valid_index()] = 0
production_MW_all.NC_obs[production_MW_all.index
                         < production_MW_all.NC_obs.first_valid_index()] = 0
# remove lines with missing values
production_MW_all = production_MW_all.dropna(axis=0)

# sum up power generation for south africa
production_MW_ZAF = pd.DataFrame({'observed': production_MW_all.EC_obs
                                  + production_MW_all.WC_obs
                                  + production_MW_all.NC_obs,
                                  'simulated': production_MW_all.EC_sim
                                  + production_MW_all.WC_sim
                                  + production_MW_all.NC_sim,
                                  'simulated_GWA': production_MW_all.EC_simGWA
                                  + production_MW_all.WC_simGWA
                                  + production_MW_all.NC_simGWA})

# drop lines before recorded generation
production_MW_ZAF = production_MW_ZAF[production_MW_ZAF.observed > 0]

# get dates as numbers
date = (production_MW_ZAF.index.year * 10000
        + production_MW_ZAF.index.month * 100
        + production_MW_ZAF.index.day)
# sum up daily and calculate from MW to GW
production_GW_ZAF_daily = production_MW_ZAF.groupby(date).sum() / 1000

# save as csv
production_GW_ZAF_daily.to_csv(results_path + "/South_Africa_daily.csv")


# -------------------------------------------------------------
# Statistical analysis for South Africa
# -------------------------------------------------------------

corZAF = [np.corrcoef(production_GW_ZAF_daily.observed,
                      production_GW_ZAF_daily.simulated)[0, 1],
          np.corrcoef(production_GW_ZAF_daily.observed,
                      production_GW_ZAF_daily.simulated_GWA)[0, 1]]
rmseZAF = [np.sqrt(np.mean((np.array(production_GW_ZAF_daily.observed)
                            - np.array(production_GW_ZAF_daily.simulated)) ** 2)),
           np.sqrt(np.mean((np.array(production_GW_ZAF_daily.observed)
                            - np.array(production_GW_ZAF_daily.simulated_GWA)) ** 2))]
mbeZAF = [(np.array(production_GW_ZAF_daily.observed)
           - np.array(production_GW_ZAF_daily.simulated)).mean(),
          (np.array(production_GW_ZAF_daily.observed)
           - np.array(production_GW_ZAF_daily.simulated_GWA)).mean()]
meanZAF = [production_GW_ZAF_daily.simulated.mean(),
           production_GW_ZAF_daily.simulated_GWA.mean(),
           production_GW_ZAF_daily.observed.mean()]

# summarise statistics
statistics_ZAF = pd.DataFrame({'South_Africa': [corZAF[0], rmseZAF[0], mbeZAF[0], meanZAF[0]],
                               'South_Africa_GWA': [corZAF[1], rmseZAF[1], mbeZAF[1], meanZAF[1]],
                               'South_Africa_observed': [None] * 3 + [meanZAF[2]]},
                              index=['Correlation', 'RMSE_GWh', 'MBE_GWh', 'Mean_GWh'])

# save results
statistics_ZAF.round(2).to_csv(results_path + "/statistics_South_Africa_daily_abs.csv", sep=";")

# plot results
boxplot3 = [list(production_GW_ZAF_daily.observed),
            list(production_GW_ZAF_daily.simulated),
            list(production_GW_ZAF_daily.simulated_GWA)]

fig3 = plt.figure(1, figsize=(9, 6))
# Create an axes instance
ax3 = fig3.add_subplot(111)
# Create the boxplot
bp3 = ax3.boxplot(boxplot3)
ax3.set_xticklabels(['observed', 'simulated', 'simulated_GWA'])


# -------------------------------------------------------------
# get mean installed capacities
# -------------------------------------------------------------

# get capacities
cap = pd.DataFrame({'capacity': list(windparks.Capacity)}, index=t)
# sort by date
cap = cap.sort_index()
# sum per date
cap = cap.groupby(cap.index).sum()
# calculate cumulative capacities
cap['capsum'] = cap.capacity.cumsum()
# get total capacity at production start
firstcap = cap[cap.index < production_MW_ZAF.index[0]].capsum[-1]
# cut capacities after starting date of production
cap_cut = cap[(cap.index >= production_MW_ZAF.index[0])
              & (cap.index <= production_MW_ZAF.index[-1])]

# plot installed capacity
plt.plot(cap_cut.capsum)

# define start and end date
startdate = np.datetime64(str(production_GW_ZAF_daily.index[0])[:4]
                          + '-' + str(production_GW_ZAF_daily.index[0])[4:6]
                          + '-' + str(production_GW_ZAF_daily.index[0])[6:])
enddate = np.datetime64(str(production_GW_ZAF_daily.index[-1])[:4]
                        + '-' + str(production_GW_ZAF_daily.index[-1])[4:6]
                        + '-' + str(production_GW_ZAF_daily.index[-1])[6:])

# prepare dataframe for fillin in capacities at commissioning dates
installed_capacity = pd.DataFrame({'capacity': np.nan}, index=np.arange(startdate, enddate))
# fill in capacities at commissioning dates
installed_capacity.capacity = installed_capacity.index.map(cap_cut.capsum)
# add capacity at start
installed_capacity.loc[installed_capacity.index[0], 'capacity'] = firstcap
# fill in gaps with last capacity
installed_capacity = installed_capacity.fillna(method='ffill')


# -------------------------------------------------------------
# use mean installed capacities to calculate relative results
# -------------------------------------------------------------
statistics_ZAF_rel = statistics_ZAF.copy(deep=True)
# factor for relativising
# division by 1000 for MW -> GW
# * 24 because day has 24 hours
rel_factor = installed_capacity.capacity.mean() / 1000 * 24
# * 100 to get result in %
statistics_ZAF_rel[1:2] = statistics_ZAF_rel[1:2] / rel_factor * 100


# save relative results
statistics_ZAF_rel.round(2).to_csv(results_path + "/statistics_South_Africa_daily_rel.csv",
                                   sep=";")
