import xarray as xr
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd

ryberg_path = "/data/users/kgruber/other-data"

def windpower_simulation_era5(winduvh,hubheight,capacity,specific_pow,lons,lats,commissioning,GWA=[]):
	'''
	function for simulating wind power generation
	
	input parameters:
	
	wind		wind speeds dataset in u and v direction in two different heights
	hubheight	vector of hub height of different turbines
	capacity	installed capacity of turbines
	lons, lats	locations of wind power plants
	commissioning	commissioning date of wind power plants
	GWA			empty list by default, if wind speed correction with GWA desired provide GWA
	'''
	
	# interpolate wind to locations of turbines
	wind = winduvh.interp(coords={"longitude":xr.DataArray(lons,dims='location'),
								  "latitude":xr.DataArray(lats,dims='location')},method="nearest")
	
	# wind in 10m height calculated from u and v component
	windh10 = (wind.v10**2+wind.u10**2)**0.5
	# wind in 100m height calculated from u and v component
	windh100 = (wind.v100**2+wind.u100**2)**0.5
	# calculate alpha friction coefficient
	alpha = (xr.ufuncs.log(windh100/windh10)/np.log(100/10))
	# calculate wind at hubheight using alpha
	windhh = windh10 * (hubheight/10)**alpha
	
	# apply GWA bias correction
	if(len(GWA)>0):
		# interpolate to turbine locations
		GWA_locations = GWA.interp(coords={"x":xr.DataArray(lons,dims='location'),
										   "y":xr.DataArray(lats,dims='location')},
								   method="nearest")
		# calculate correction factor
		cf_GWA = GWA_locations.sel(band=1)/windh100.mean('time')
		# apply correction factor
		windhh = windhh * cf_GWA
	
	# replace wind speeds higher than 25 m/s with 0, because cutout windspeed
	windhh = windhh.where(windhh<=25,0)
	
	# Ryberg power curve model
	RybCoeff = pd.read_csv(ryberg_path+"/ryberg_coeff.csv")
	A = xr.DataArray(RybCoeff.A, dims = 'CF')
	B = xr.DataArray(RybCoeff.B, dims = 'CF')
	ignore_dims = [[], [], ['CF'], ['CF']]
	wp1 = xr.apply_ufunc(power_curve,windhh,specific_pow,A,B,
                         input_core_dims = ignore_dims,
                         dask = 'parallelized',
                         output_dtypes = [np.float64],
                         vectorize = True)

	# multiply with installed capacity
	wp2 = capacity*wp1
	
	# make wind power generation start at commissioning date
	if(len(GWA)>0):
		wp3 =  wp2.where(wp2.time >= xr.DataArray(commissioning,dims={'location':range(len(commissioning))}), 0).compute().drop(['x','y','band'])
	else:
		wp3 =  wp2.where(wp2.time >= xr.DataArray(commissioning,dims={'location':range(len(commissioning))}), 0).compute()
	
	return(wp3)


# function for calculating power output with Ryberg model
# wind: wind speeds over time at one location
# spec_pow: specific power (W/m^2) for that turbine/park/location
# A, B: coefficients of Ryberg model
def power_curve(wind,spec_pow,A,B):
    v = [0] + list(np.exp(A+B*np.log(spec_pow))) + [25]
    CF = [0] + list(range(101)) + [100]
    PC = interp1d(v, CF)
    return(PC(wind))

	
def windpower_simulation_merra2(winduvh,hubheight,capacity,specific_pow,lons,lats,commissioning,GWA=[]):
	'''
	function for simulating wind power generation
	
	input parameters:
	
	wind		wind speeds dataset in u and v direction in two different heights
	hubheight	vector of hub height of different turbines
	capacity	installed capacity of turbines
	lons, lats	locations of wind power plants
	commissioning	commissioning date of wind power plants
	GWA			empty list by default, if wind speed correction with GWA desired provide GWA
	'''
	
	# interpolate wind to locations of turbines
	wind = winduvh.interp(coords={"lon":xr.DataArray(lons,dims='location'),
								  "lat":xr.DataArray(lats,dims='location')},method="nearest")
	
	# wind in 10m height calculated from u and v component
	windh10 = (wind.V10M**2+wind.U10M**2)**0.5
	# wind in 100m height calculated from u and v component
	windh50 = (wind.V50M**2+wind.U50M**2)**0.5
	# calculate alpha friction coefficient
	alpha = (xr.ufuncs.log(windh50/windh10)/np.log(50/(10+wind.DISPH)))
	# calculate wind at hubheight using alpha
	windhh = windh50 * (hubheight/50)**alpha
	
	# apply GWA bias correction
	if(len(GWA)>0):
		# interpolate to turbine locations
		GWA_locations = GWA.interp(coords={"x":xr.DataArray(lons,dims='location'),
										   "y":xr.DataArray(lats,dims='location')},
								   method="nearest")
		# calculate correction factor
		cf_GWA = GWA_locations.sel(band=1)/windh50.mean('time')
		# apply correction factor
		windhh = windhh * cf_GWA
	
	# replace wind speeds higher than 25 m/s with 25, because maximum of power curve
	windhh = windhh.where(windhh<=25,25)
	
	# Ryberg power curve model
	RybCoeff = pd.read_csv(ryberg_path+"/ryberg_coeff.csv")
	A = xr.DataArray(RybCoeff.A, dims = 'CF')
	B = xr.DataArray(RybCoeff.B, dims = 'CF')
	ignore_dims = [[], [], ['CF'], ['CF']]
	wp1 = xr.apply_ufunc(power_curve,windhh,specific_pow,A,B,
                         input_core_dims = ignore_dims,
                         dask = 'parallelized',
                         output_dtypes = [np.float64],
                         vectorize = True)

	# multiply with installed capacity
	wp2 = capacity*wp1
	
	
	# make wind power generation start at commissioning date
	if(len(GWA)>0):
		wp3 =  wp2.where(wp2.time >= xr.DataArray(commissioning,dims={'location':range(len(commissioning))}), 0).compute().drop(['x','y','band'])
	else:
		wp3 =  wp2.where(wp2.time >= xr.DataArray(commissioning,dims={'location':range(len(commissioning))}), 0).compute()
	
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