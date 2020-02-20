import glob
import pandas as pd
import numpy as np
import os
import geopandas

from paths_bra import *

turb_path = bra_path + '/aerogeradores'

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

# merge BIG commissioning dates into turbine data
com_anl = pd.Series(BIG_anl['Data Operação'].values,index=BIG_anl.CEG2)
wp_anl['commissioning'] = wp_anl.CEG.map(com_anl)

# not all codes match - match remaining turbines by wind park names
com_anlN = pd.Series(BIG_anl['Data Operação'].values,index=BIG_anl.Usina)
wp_anl.loc[pd.isnull(wp_anl.commissioning),'commissioning'] = wp_anl.name[pd.isnull(wp_anl.commissioning)].map(com_anlN)

# convert commissioning to datetime format
wp_anl['commissioning'] = wp_anl.commissioning.replace('-',np.nan)
wp_anl.commissioning = pd.to_datetime(wp_anl.commissioning)
# make capacities numeric
wp_anl['cap'] = wp_anl.cap.astype(np.float)*1000

wp_anl.to_csv(bra_path + '/windturbines_BRA_complete.csv')