#!/usr/bin/bash

# download data
python download_era5_ZAF.py
echo "era5 downloaded"
python download_merra2_ZAF.py
echo "merra2 downloaded"

# call scripts for preparing reanalysis data
python prepare_ZAF_ERA5.py
echo "era5 prepared"
python prepare_ZAF_MERRA2.py
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

# cut out GWA2
python cut_GWA2.py

DATASETS='MERRA ERA'
for DS in $DATASETS
do
    echo "simulate " $DS
    python WP_simulation_ZAF_$DS.py
    echo "simulate " $DS "GWA"
    python WP_simulation_ZAF_${DS}_GWA.py -GWA $GWA
done

# analyse results
python WP_simulation_ZAF_Analysis.py
# calculate grid size
python number_reanalysis_grid_points.py
