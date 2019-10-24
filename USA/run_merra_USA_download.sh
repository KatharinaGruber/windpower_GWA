#!/usr/bin/bash
N=${1:-2}

MONTHS='01 02 03 04 05 06 07 08 09 10 11 12'

for year in {2000..2018}
do
   for month in $MONTHS
   do
      sem --ungroup -j $N --ungroup python download_merra2_USA.py -date $year$month
   done
done
sem --wait