import argparse
import geopandas
import glob
import numpy as np
import os
import pandas as pd
from functools import reduce
import xarray as xr

from paths_bra import *

turb_path = bra_path + '/aerogeradores'

parser = argparse.ArgumentParser(description='Insert optionally GWA')
parser.add_argument('-GWA')
args = parser.parse_args()
if(args.GWA == None):
    GWA = "3"
else:
    GWA = args.GWA

if bra_path + "/turbine_data.csv" not in glob.glob(bra_path + '/*.csv'):
	# special characters dictionary
	specials = {'Ã§':'ç',
		     'Ã£':'ã',
		     'Ã\xad':'í',
		     'Ã¡':'á',
		     'Ã³':'ó',
		     'Ã©':'é',
		     'Ã¢':'â',
		     'Ãº':'ú',
		     'Ã´':'ô',
		     'Ãª':'ê',
		     'Ã\x81':'Á'}


	# prepare turbine data: generate csv from shapefile
	turb_file = turb_path + '/Aerogeradores.csv'
	if turb_file not in glob.glob(turb_path + '/*'):
	    gdf = geopandas.read_file(turb_path + '/Aerogeradores.shp').to_crs({'init': 'epsg:4326'})
	    pd.DataFrame({'cap': gdf.POT_MW,
		          'height': gdf.ALT_TORRE,
		          'diam': gdf.DIAM_ROTOR,
		          'name': gdf.NOME_EOL,
		          'CEG':gdf.CEG,
		          'operating': gdf.Operacao,
		          'lon': gdf.geometry.x,
		          'lat': gdf.geometry.y}).to_csv(turb_file)

	# read turbine data and select only operating turbines
	wp_anl = pd.read_csv(turb_file,index_col=0)
	wp_anl = wp_anl[wp_anl.operating=='SIM'].drop('operating',axis=1)
	for i in specials:
	    wp_anl.name = wp_anl.name.str.replace(i,specials.get(i)) # correct special characters
	wp_anl['name'] = wp_anl.name.str.replace('Umburanas 0','Umburanas ').str.replace('Zeus','ZEUS') # adapt slightly different names

	# load BIG data and prepare CEG code to match format in turbine data
	BIG_anl = pd.read_csv(bra_path + '/BIG_windparks.csv',header=2)[:-7] # remove last 7 rows with additional information
	BIG_anl['CEG2'] = BIG_anl.CEG.str.replace(r'[.]', '',n=3).str.replace(r'[.]','-') # adapt codes to be the same format
	# add state
	BIG_anl['state'] = BIG_anl.CEG.str.split('.',expand=True)[2]

	# merge BIG commissioning dates and state into turbine data
	com_anl = pd.DataFrame({'com': BIG_anl['Data Operação'].values,
		                 'state': BIG_anl.state.values}, index=BIG_anl.CEG2)
	wp_anl['commissioning'] = wp_anl.CEG.map(com_anl.com)
	wp_anl['state'] = wp_anl.CEG.map(com_anl.state)

	# not all codes match - match remaining turbines by wind park names
	com_anlN = pd.DataFrame({'com': BIG_anl['Data Operação'].values,
		                 'state': BIG_anl.state.values}, index=BIG_anl.Usina)
	wp_anl.loc[pd.isnull(wp_anl.commissioning),'commissioning'] = wp_anl.name[pd.isnull(wp_anl.commissioning)].map(com_anlN.com)
	wp_anl.loc[pd.isnull(wp_anl.state),'state'] = wp_anl.name[pd.isnull(wp_anl.state)].map(com_anlN.state)

	# convert commissioning to datetime format
	wp_anl['commissioning'] = wp_anl.commissioning.replace('-',np.nan)
	wp_anl.commissioning = pd.to_datetime(wp_anl.commissioning)
	# make capacities numeric
	wp_anl['cap'] = wp_anl.cap.astype(np.float)*1000

	# two wind parks have 0 hub heights and diameters
	# -> replace them by diameter of similar windturbines (same capacity and commissioning year)
	# diameter: 100 m, height: 95 m
	# there also is one windpark with little distance of rotor blades to ground (Esperança)
	# - is this a problem?
	wp_anl['diam'][wp_anl.diam==0] = 100
	wp_anl['height'][wp_anl.height==0] = 95

	# calculate specific power
	wp_anl['sp'] = wp_anl.cap*1000/(wp_anl.diam**2/4*np.pi)


	#wp_anl.to_csv(bra_path + '/windturbines_BRA_complete.csv')
	wp_anl.to_csv(bra_path + '/turbine_data.csv')




if bra_path + '/turbine_data_mer_gwa' + GWA + '.csv' not in glob.glob(bra_path + '/*.csv'):

	# prepare turbine data labels
	turbine_data = pd.read_csv(bra_path + "/turbine_data.csv", index_col = [0])
	turbine_data['ind'] = range(len(turbine_data))

	# MERRA2
	data_mer = xr.open_dataset(mer_path+"/merra2_wind_BRA_200601.nc")
	# Create dataframe with sequence the size of MERRA-2 grid to find out which turbines interpolate to the same point
	in_seq = xr.Dataset({'x':(['lat','lon'],
		                  np.array(range(data_mer.DISPH.isel(time=0).values.size)).reshape(data_mer.DISPH.isel(time=0).values.shape))},
		             coords = {'lat':data_mer.lat.values,
		                       'lon':data_mer.lon.values})
	# interpolate to indices
	ip = in_seq.interp(coords={"lon":xr.DataArray(turbine_data.lon,dims='location'),
		                   "lat":xr.DataArray(turbine_data.lat,dims='location')},method="nearest").to_dataframe()
	# Load GWA data and extract values at locations to find unique locations
	if GWA == "3":
		GWA_BRA = xr.open_rasterio(bra_path+'/GWA/GWA3_BRA50m.tif')
	else:
		GWA_BRA = xr.open_rasterio(bra_path+'/GWA/GWA_BRA50m.tif')
	# interpolate GWA to all locations and compare if same locations from rea dataset match same locations from GWA
	GWA_ip = GWA_BRA.sel(band=1).interp(coords={"x":xr.DataArray(turbine_data.lon,dims='location'),
												"y":xr.DataArray(turbine_data.lat,dims='location')},
										method="nearest").to_dataframe(name='GWA')
	# create labels for accumulating data to determine which configuartions are unique label consists of:
	# - state
	# - interpolation point
	# - commissioning date
	# - specific power
	# - hub height
	lbl_mer = reduce(np.core.defchararray.add,[turbine_data.state.values.astype('U'),'_',
											   np.array(ip.x).astype('U'),'_',
											   turbine_data.commissioning.astype('U'),"_sp",
											   turbine_data.sp.astype('U'),"_hh",
											   turbine_data.height.astype('U')])
	lbl_mer_gwa = reduce(np.core.defchararray.add,[turbine_data.state.values.astype('U'),'_',
												   np.array(ip.x).astype('U'),'_',
												   turbine_data.commissioning.astype('U'),"_sp",
												   turbine_data.sp.astype('U'),"_hh",
												   turbine_data.height.astype('U'), "_gwa",
												   GWA_ip.GWA.astype('U')])

	state_ag = turbine_data.state.groupby(lbl_mer).first().values
	cap_ag = turbine_data.cap.groupby(lbl_mer).sum().values
	ind_ag = turbine_data.ind.groupby(lbl_mer).min().values
	comm_ag = turbine_data.commissioning.values[ind_ag]
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
	state_ag = turbine_data.state.groupby(lbl_mer_gwa).first().values
	cap_ag = turbine_data.cap.groupby(lbl_mer_gwa).sum().values
	ind_ag = turbine_data.ind.groupby(lbl_mer_gwa).min().values
	comm_ag = turbine_data.commissioning.values[ind_ag]
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
		                             

	# ERA5
	data_era = xr.open_dataset(era_path+"/era5_wind_BRA_200601.nc")
	# Create dataframe with sequence the size of MERRA-2 grid to find out which turbines interpolate to the same point
	in_seq = xr.Dataset({'x':(['lat','lon'],
		                  np.array(range(data_era.u10.isel(time=0).values.size)).reshape(data_era.u10.isel(time=0).values.shape))},
		             coords = {'lat':data_era.latitude.values,
		                       'lon':data_era.longitude.values})
	# interpolate to indices
	ip = in_seq.interp(coords={"lon":xr.DataArray(turbine_data.lon,dims='location'),
		                   "lat":xr.DataArray(turbine_data.lat,dims='location')},method="nearest").to_dataframe()
	# Load GWA data and extract values at locations to find unique locations	
	if GWA == "3":    
		GWA_BRA = xr.open_rasterio(bra_path+'/GWA/GWA3_BRA100m.tif')
	else:
		GWA_BRA = xr.open_rasterio(bra_path+'/GWA/GWA_BRA100m.tif')
	# interpolate GWA to all locations and compare if same locations from rea dataset match same locations from GWA
	GWA_ip = GWA_BRA.sel(band=1).interp(coords={"x":xr.DataArray(turbine_data.lon,dims='location'),
		                                    "y":xr.DataArray(turbine_data.lat,dims='location')},
		                            method="nearest").to_dataframe(name='GWA')
	# create labels for accumulating data to determine which configuartions are unique label consists of:
	# - state
	# - interpolation point
	# - commissioning date
	# - specific power
	# - hub height
	lbl_era = reduce(np.core.defchararray.add,[turbine_data.state.values.astype('U'),'_',
		                                   turbine_data.state.astype('U'),'_',
		                                   np.array(ip.x).astype('U'),'_',
		                                   turbine_data.commissioning.astype('U'),"_sp",
		                                   turbine_data.sp.astype('U'),"_hh",
		                                   turbine_data.height.astype('U')])
	lbl_era_gwa = reduce(np.core.defchararray.add,[turbine_data.state.values.astype('U'),'_',
		                                       turbine_data.state.astype('U'),'_',
		                                       np.array(ip.x).astype('U'),'_',
		                                       turbine_data.commissioning.astype('U'),"_sp",
		                                       turbine_data.sp.astype('U'),"_hh",
		                                       turbine_data.height.astype('U'), "_gwa",
		                                       GWA_ip.GWA.astype('U')])
	comm_state = turbine_data.state.groupby(lbl_era).first().values
	cap_ag = turbine_data.cap.groupby(lbl_era).sum().values
	ind_ag = turbine_data.ind.groupby(lbl_era).min().values
	comm_ag = turbine_data.commissioning.values[ind_ag]
	lon_ag = turbine_data.lon.values[ind_ag]
	lat_ag = turbine_data.lat.values[ind_ag]
	sp_ag = turbine_data.sp.values[ind_ag]
	hh_ag = turbine_data.height.values[ind_ag]
	turbine_data_era = pd.DataFrame({'state':comm_state,
		                         'commissioning':comm_ag,
		                         'capacity':cap_ag,
		                         'lon':lon_ag,
		                         'lat':lat_ag,
		                         'sp':sp_ag,
		                         'height':hh_ag})
	# with GWA
	comm_state = turbine_data.state.groupby(lbl_era_gwa).first().values
	cap_ag = turbine_data.cap.groupby(lbl_era_gwa).sum().values
	ind_ag = turbine_data.ind.groupby(lbl_era_gwa).min().values
	comm_ag = turbine_data.commissioning.values[ind_ag]
	lon_ag = turbine_data.lon.values[ind_ag]
	lat_ag = turbine_data.lat.values[ind_ag]
	sp_ag = turbine_data.sp.values[ind_ag]
	hh_ag = turbine_data.height.values[ind_ag]
	turbine_data_era_gwa = pd.DataFrame({'state':comm_state,
		                             'commissioning':comm_ag,
		                             'capacity':cap_ag,
		                             'lon':lon_ag,
		                             'lat':lat_ag,
		                             'sp':sp_ag,
		                             'height':hh_ag})
		                             
	# save to files
	turbine_data_mer.to_csv(bra_path + '/turbine_data_mer.csv')
	turbine_data_mer_gwa.to_csv(bra_path + '/turbine_data_mer_gwa' + GWA + '.csv')
	turbine_data_era.to_csv(bra_path + '/turbine_data_era.csv')
	turbine_data_era_gwa.to_csv(bra_path + '/turbine_data_era_gwa' + GWA + '.csv')

	# merge labels and save
	labels = pd.DataFrame({'lbl_mer':lbl_mer,
		               'lbl_mer_gwa':lbl_mer_gwa,
		               'lbl_era':lbl_era,
		               'lbl_era_gwa':lbl_era_gwa})
	labels.to_csv(bra_path + '/labels_turbine_data_gwa' + GWA + '.csv')