#!/bin/bash

## this script creates folders for separate simulations in a given temperature range

if [ "$#" -ne 5 ]; then
    echo "Usage: <script name> Tmin Tmax dT beta_SU2 r_U1"
	exit 1
fi

Tmin=$1
Tmax=$2
dT=$3
beta=$4
ru1=$5



for T in $(seq $Tmin $dT $Tmax);
do

	DIRNAME="T-$T"

	if [ ! -d "$DIRNAME" ]; then
  		mkdir $DIRNAME
		write_config_singlet.sh ../input_couplings $T $beta $ru1
		cp config $DIRNAME
		mv params_lattice.dat $DIRNAME
		mv params_continuum.dat $DIRNAME
		#sed -i "s/--job-name=.*/--job-name=T$T/g" submit.job
		cp submit.job $DIRNAME
	else
		echo 'Error: directory' $DIRNAME 'exists.'
	fi


done

