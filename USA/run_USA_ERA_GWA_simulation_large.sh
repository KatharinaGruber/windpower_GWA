#!/usr/bin/bash
N=${1:-2}
GWA=$2
STATES='CA CO IA IL KS OK TX'
for state in $STATES
do
   sem --ungroup -j $N --ungroup python WP_simulation_USA_ERA_GWA_large.py -state $state -GWA $GWA
done
sem --wait
