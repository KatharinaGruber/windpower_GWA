#!/usr/bin/bash
N=${1:-2}
GWA=$2
STATES='BA CE MA PB PE PI PR RJ RN RS SC SE'
for state in $STATES
do
   sem --ungroup -j $N --ungroup python WP_simulation_BRA_ERA_GWA.py -state $state -GWA $GWA
done
sem --wait
