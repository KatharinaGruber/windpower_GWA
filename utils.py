import xarray as xr
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd

from paths import ryberg_path

def windpower_simulation_era5(windh100,alpha,hubheight,capacity,lons,lats,commissioning,GWA=[]):
    '''
    function for simulating wind power generation
    
    input parameters:
    
    windh100	wind speeds dataset with effective wind speeds in one height (100m)
    alpha		alpha friction coefficient calculate from wind speeds in two heights
    hubheight	vector of hub height of different turbines
    capacity	installed capacity of turbines
    lons, lats	locations of wind power plants
    commissioning	commissioning date of wind power plants
    GWA			empty list by default, if wind speed correction with GWA desired provide GWA
    '''
    
    # interpolate wind to locations of turbines
    wind = windh100.interp(coords={"longitude":xr.DataArray(lons,dims='location'),
                                   "latitude":xr.DataArray(lats,dims='location')},method="nearest").compute()
    # interpolate alpha to locations of turbines
    alphai = alpha.interp(coords={"longitude":xr.DataArray(lons,dims='location'),
                                  "latitude":xr.DataArray(lats,dims='location')},method="nearest").compute()
    # calculate wind at hubheight using alpha
    windhh = (wind * (hubheight/100)**alphai).compute()
    
    # apply GWA bias correction
    if(len(GWA)>0):
        # interpolate to turbine locations
        GWA_locations = GWA.interp(coords={"x":xr.DataArray(lons,dims='location'),
                                           "y":xr.DataArray(lats,dims='location')},
                                   method="nearest").compute()
        # calculate correction factor
        cf_GWA = (GWA_locations.sel(band=1)/wind.mean('time')).compute()
        # apply correction factor
        windhhg = (windhh * cf_GWA).compute()
        # replace wind speeds higher than 25 m/s with 0, because cutout windspeed
        windhhg = windhhg.where(windhhg<=25,0)
    else:
        windhhg = windhh.where(windhh<=25,0)
    
    # Enercon E-82 power curve
	wind_speeds = (np.arange(0, 26, step=1.0))
	generation_kw = [0.0, 0.000000000001, 3.0, 25.0, 82.0, 175.0, 321.0, 532.0, 815.0, 1180.0, 1580.0, 1810.0, 1980.0] + 13 * [2050.0]
	
	power_curve = interp1d(wind_speeds, generation_kw)
	
	wp1 = xr.apply_ufunc(power_curve, windhh,
                     dask='parallelized',
                     output_dtypes=[np.float64])
	# fetch installed capacity and divide by 2000 to make factor for capacity of Enercon E-82
	cap = list(capacity/2000.0)
	# multiply with installed capacity
	wp2 = cap*wp1
    
    # make wind power generation start at commissioning date
    if(len(GWA)>0):
        wp3 =  wp2.where(wp2.time >= xr.DataArray(commissioning,coords={'location':range(len(commissioning))},dims='location'), 0).compute().drop('band')
    else:
        wp3 =  wp2.where(wp2.time >= xr.DataArray(commissioning,coords={'location':range(len(commissioning))},dims='location'), 0).compute()
    
    return(wp3)




    
def windpower_simulation_merra2(windh50,alpha,hubheight,capacity,lons,lats,commissioning,GWA=[]):
    '''
    function for simulating wind power generation
    
    input parameters:
    
    windh50		wind speeds dataset with effective wind speeds in one height (50m)
    alpha		alpha friction coefficient calculate from wind speeds in two heights
    hubheight	vector of hub height of different turbines
    capacity	installed capacity of turbines
    lons, lats	locations of wind power plants
    commissioning	commissioning date of wind power plants
    GWA			empty list by default, if wind speed correction with GWA desired provide GWA
    '''
    
    # interpolate wind to locations of turbines
    wind = windh50.interp(coords={"lon":xr.DataArray(lons,dims='location'),
                                  "lat":xr.DataArray(lats,dims='location')},method="nearest").compute()
    
    alphai = alpha.interp(coords={"lon":xr.DataArray(lons,dims='location'),
                                  "lat":xr.DataArray(lats,dims='location')},method="nearest").compute()
    # calculate wind at hubheight using alpha
    windhh = (wind * (hubheight/50)**alphai).compute()
    
    # apply GWA bias correction
    if(len(GWA)>0):
        # interpolate to turbine locations
        GWA_locations = GWA.interp(coords={"x":xr.DataArray(lons,dims='location'),
                                           "y":xr.DataArray(lats,dims='location')},
                                   method="nearest").compute()
        # calculate correction factor
        cf_GWA = (GWA_locations.sel(band=1)/wind.mean('time')).compute()
        # apply correction factor
        windhhg = windhh * cf_GWA
        # replace wind speeds higher than 25 m/s with 25, because maximum of power curve
        windhhg = windhhg.where(windhhg<=25,0)
    else:
        windhhg = windhh.where(windhh<=25,0)
    
    # Enercon E-82 power curve
	wind_speeds = (np.arange(0, 26, step=1.0))
	generation_kw = [0.0, 0.000000000001, 3.0, 25.0, 82.0, 175.0, 321.0, 532.0, 815.0, 1180.0, 1580.0, 1810.0, 1980.0] + 13 * [2050.0]
	
	power_curve = interp1d(wind_speeds, generation_kw)
	
	wp1 = xr.apply_ufunc(power_curve, windhh,
                     dask='parallelized',
                     output_dtypes=[np.float64])
	# fetch installed capacity and divide by 2000 to make factor for capacity of Enercon E-82
	cap = list(capacity/2000.0)
	# multiply with installed capacity
	wp2 = cap*wp1
    
    
    # make wind power generation start at commissioning date
    if(len(GWA)>0):
        wp3 =  wp2.where(wp2.time >= xr.DataArray(commissioning,coords={'location':range(len(commissioning))},dims='location'), 0).compute().drop('band')
    else:
        wp3 =  wp2.where(wp2.time >= xr.DataArray(commissioning,coords={'location':range(len(commissioning))},dims='location'), 0).compute()
    
    return(wp3)
    
    
    
def stats(sim,obs,rd=True):
    '''
    function for calculating statistics for analysis of simulated wind power generation time series
        
    input parameters:
    
    sim		simulated windpower time series, pandas Series
    obs		observed windpower time series, pandas Series
    rd		round results to two digits
    
    sim and obs need to have the same length and index
    
    '''
    cor = np.corrcoef(sim,obs)[0,1]
    rmse = (((sim-obs)**2).mean())**0.5
    mbe = (sim-obs).mean()
    mean = sim.mean()
    
    if rd:
        return([round(i,2) for i in [cor,rmse,mbe,mean]])
    else:
        return([cor,rmse,mbe,mean])
        
        
        
          
    
    
# functions for handling larger states (with more locations)
def windpower_simulation_era5_large(windh100,alpha,hubheight,capacity,lons,lats,commissioning,GWA=[]):
    # interpolate wind to locations of turbines
    wind = windh100.interp(coords={"longitude":xr.DataArray(lons,dims='location'),
                                   "latitude":xr.DataArray(lats,dims='location')},method="nearest")
    # interpolate alpha to locations of turbines
    alphai = alpha.interp(coords={"longitude":xr.DataArray(lons,dims='location'),
                                  "latitude":xr.DataArray(lats,dims='location')},method="nearest")
    # calculate wind at hubheight using alpha
    windhh = (wind * (hubheight/100)**alphai)

    # apply GWA bias correction
    if(len(GWA)>0):
        # interpolate to turbine locations
        GWA_locations = GWA.interp(coords={"x":xr.DataArray(lons,dims='location'),
                                           "y":xr.DataArray(lats,dims='location')},
                                   method="nearest")
        # calculate correction factor
        cf_GWA = (GWA_locations.sel(band=1)/wind.mean('time'))
        # apply correction factor
        windhhg = (windhh * cf_GWA)
        # replace wind speeds higher than 25 m/s with 0, because cutout windspeed
        windhhg = windhhg.where(windhhg<=25,0)
    else:
        windhhg = windhh.where(windhh<=25,0)

    # Enercon E-82 power curve
	wind_speeds = (np.arange(0, 26, step=1.0))
	generation_kw = [0.0, 0.000000000001, 3.0, 25.0, 82.0, 175.0, 321.0, 532.0, 815.0, 1180.0, 1580.0, 1810.0, 1980.0] + 13 * [2050.0]
	
	power_curve = interp1d(wind_speeds, generation_kw)
	
	wp1 = xr.apply_ufunc(power_curve, windhh,
                     dask='parallelized',
                     output_dtypes=[np.float64])
	# fetch installed capacity and divide by 2000 to make factor for capacity of Enercon E-82
	cap = list(capacity/2000.0)
	# multiply with installed capacity
	wp2 = cap*wp1

    # make wind power generation start at commissioning date
    if(len(GWA)>0):
        wp3 =  wp2.where(wp2.time >= xr.DataArray(commissioning,coords={'location':range(len(commissioning))},dims='location'), 0).compute().drop('band')
    else:
        wp3 =  wp2.where(wp2.time >= xr.DataArray(commissioning,coords={'location':range(len(commissioning))},dims='location'), 0).compute()

    return(wp3)
    
def windpower_simulation_merra2_large(windh50,alpha,hubheight,capacity,lons,lats,commissioning,GWA=[]):
    # interpolate wind to locations of turbines
    wind = windh50.interp(coords={"lon":xr.DataArray(lons,dims='location'),
                                  "lat":xr.DataArray(lats,dims='location')},method="nearest")
    # interpolate alpha to locations of turbines
    alphai = alpha.interp(coords={"lon":xr.DataArray(lons,dims='location'),
                                  "lat":xr.DataArray(lats,dims='location')},method="nearest")
    # calculate wind at hubheight using alpha
    windhh = (wind * (hubheight/50)**alphai)

    # apply GWA bias correction
    if(len(GWA)>0):
        # interpolate to turbine locations
        GWA_locations = GWA.interp(coords={"x":xr.DataArray(lons,dims='location'),
                                           "y":xr.DataArray(lats,dims='location')},
                                   method="nearest")
        # calculate correction factor
        cf_GWA = (GWA_locations.sel(band=1)/wind.mean('time'))
        # apply correction factor
        windhhg = (windhh * cf_GWA)
        # replace wind speeds higher than 25 m/s with 0, because cutout windspeed
        windhhg = windhhg.where(windhhg<=25,0)
    else:
        windhhg = windhh.where(windhh<=25,0)

    # Enercon E-82 power curve
	wind_speeds = (np.arange(0, 26, step=1.0))
	generation_kw = [0.0, 0.000000000001, 3.0, 25.0, 82.0, 175.0, 321.0, 532.0, 815.0, 1180.0, 1580.0, 1810.0, 1980.0] + 13 * [2050.0]
	
	power_curve = interp1d(wind_speeds, generation_kw)
	
	wp1 = xr.apply_ufunc(power_curve, windhh,
                     dask='parallelized',
                     output_dtypes=[np.float64])
	# fetch installed capacity and divide by 2000 to make factor for capacity of Enercon E-82
	cap = list(capacity/2000.0)
	# multiply with installed capacity
	wp2 = cap*wp1

    # make wind power generation start at commissioning date
    if(len(GWA)>0):
        wp3 =  wp2.where(wp2.time >= xr.DataArray(commissioning,coords={'location':range(len(commissioning))},dims='location'), 0).compute().drop('band')
    else:
        wp3 =  wp2.where(wp2.time >= xr.DataArray(commissioning,coords={'location':range(len(commissioning))},dims='location'), 0).compute()

    return(wp3)
	
