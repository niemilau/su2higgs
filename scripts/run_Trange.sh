#!/bin/bash

## this script runs a single-node simulation in each T-folder (so run prepare_Trange.sh first)

if [ "$#" -ne 3 ]; then
    echo "Usage: <script name> Tmin Tmax dT"
	exit 1
fi

Tmin=$1
Tmax=$2
dT=$3


for T in $(seq $Tmin $dT $Tmax);
do

	DIRNAME="T-$T"

	if [ ! -d "$DIRNAME" ]; then
  		echo 'Error: No directory' $DIRNAME 
	else
		cd $DIRNAME
		condor_submit submit.job
		sleep 0.1
		cd .. 
	fi

done

