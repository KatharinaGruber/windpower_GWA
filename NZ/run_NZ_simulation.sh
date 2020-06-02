#!/usr/bin/bash

# download data
python download_era5_NZ.py
echo "era5 downloaded"
python download_merra2_NZ.py
echo "merra2 downloaded"

# download generation data
python download_generation_data.py

# call scripts for preparing reanalysis data
python prepare_NZ_ERA5.py
echo "era5 prepared"
python prepare_NZ_MERRA2.py
echo "merra2 prepared"

# run simulations without and with GWA
# determine which GWA
# if no info given -> use GWA3 as standard
if [ ! -z $1 ] 
then 
    GWA=$1
else
    GWA="3"
fi

DATASETS='MERRA ERA'
for DS in $DATASETS
do
    echo "simulate " $DS
    python WP_simulation_NZ_$DS.py
    echo "simulate " $DS "GWA"
    python WP_simulation_NZ_${DS}_GWA.py -GWA $GWA
done

#python WP_simulation_BRA_Analysis.py -GWA $GWA

