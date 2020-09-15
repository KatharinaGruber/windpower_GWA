#!/usr/bin/bash

# download data
python download_era5_BRA.py
echo "era5 downloaded"
python download_merra2_BRA.py
echo "merra2 downloaded"

# call scripts for preparing reanalysis data
python prepare_BRA_ERA5.py
echo "era5 prepared"
python prepare_BRA_MERRA2.py
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

if [ $GWA = "2" ]
then
	echo "cut out GWA2"
	python cut_GWA2.py
fi

# prepare turbine data
python prepare_BRA_turbines.py -GWA $GWA

DATASETS='MERRA ERA'
for DS in $DATASETS
do
    echo "simulate " $DS
    python WP_simulation_BRA_$DS.py
    echo "simulate " $DS "GWA"
    bash run_BRA_${DS}_GWA_simulation.sh 1 $GWA
done

python WP_simulation_BRA_Analysis.py -GWA $GWA

