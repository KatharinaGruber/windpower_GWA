#!/usr/bin/bash
N=${1:-2}
STATES='AK AR AZ CA CO CT DE FL HI IA ID IL IN KS MA MD ME MI MN MO MT NC ND NE NH NJ NM NV NY OH OK OR PA PR RI SD TN TX UT VT WA WI WV WY'
for state in $STATES
do
   sem --ungroup -j $N --ungroup python sim_onestate_old_ERAGWA.py -state $state
done
sem --wait