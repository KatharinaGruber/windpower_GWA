#!/usr/bin/bash

# download data
python download_era5_USA.py
echo "era5 downloaded"
python download_merra2_USA.py
echo "merra2 downloaded"

# call scripts for preparing reanalysis data
python prepare_USA_ERA5.py
echo "era5 prepared"
python prepare_USA_MERRA2.py
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

# prepare turbine data
python prepare_USA_turbines.py -GWA $GWA

DATASETS='MERRA ERA'
for DS in $DATASETS
do
    echo "simulate " $DS
    python WP_simulation_USA_$DS.py
    echo "simulate " $DS "GWA"
    bash run_USA_${DS}_GWA_simulation.sh 1 $GWA
    echo "simulate " $DS "GWA large"
    bash run_USA_${DS}_GWA_simulation_large.sh 1 $GWA
done

python sum_up_regions_USA.py -GWA $GWA

python WP_simulation_USA_Analysis.py -GWA $GWA

python WP_simulation_USA_Analysis2010.py -GWA $GWA
