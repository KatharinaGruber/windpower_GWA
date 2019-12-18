#!/usr/bin/bash

# download data
python download_era5_USA.py
echo "era5 downloaded"
bash bash run_merra_USA_download.sh 8
echo "merra2 downloaded"

# call scripts for preparing reanalysis data
python prepare_USA_ERA5.py
echo "era5 prepared"
python prepare_USA_MERRA2.py
echo "merra2 prepared"

# run simulations without and with GWA
DATASETS='MERRA ERA'
for DS in $DATASETS
do
    echo "simulate " $DS
	python WP_simulation_USA_$DS.py
	echo "simulate " $DS "GWA"
	bash run_USA_${DS}_GWA_simulation.sh 1
	echo "simulate " $DS "GWA large"
	bash run_USA_${DS}_GWA_simulation_large.sh 1
done

python sum_up_regions_USA.py

python WP_simulation_USA_Analysis.py
