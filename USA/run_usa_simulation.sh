#!/usr/bin/bash

# download data
#python download_era5_USA.py
#bash bash run_merra_USA_download.sh 8

# call scripts for preparing reanalysis data
#python prepare_USA_ERA5.py
#python prepare_USA_MERRA2.py

# run simulations without and with GWA
DATASETS='MERRA ERA'
for DS in $DATASETS
do
    python WP_simulation_USA_$DS.py
	bash run_USA_{$DS}_GWA_simulation.sh 1
	bash run_USA_{$DS}_GWA_simulation_large.sh 1
done

python sum_up_regions.py

python WP_simulation_USA_Analysis.py
