#!/usr/bin/bash
N=${1:-2}
STATES='CA IA KS OK TX'
for state in $STATES
do
   sem --ungroup -j $N --ungroup python WP_simulation_USA_MERRA_GWA_large.py -state $state
done
sem --wait