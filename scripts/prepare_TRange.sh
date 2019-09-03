#!/bin/bash

## this script creates folders for separate simulations in a given temperature range

if [ "$#" -ne 4 ]; then
    echo "Usage: <script name> Tmin Tmax dT betaG"
	exit 1
fi

Tmin=$1
Tmax=$2
dT=$3
beta=$4



for T in $(seq $Tmin $dT $Tmax);
do

	DIRNAME="T-$T"

	if [ ! -d "$DIRNAME" ]; then
  		mkdir $DIRNAME
		cp su2 $DIRNAME
		./write_config.sh input_couplings $T $beta
		cp config $DIRNAME
		mv params_lattice.dat $DIRNAME
		mv params_MSbar.dat $DIRNAME
	else
		echo 'Error: directory' $DIRNAME 'exists.'
	fi

	echo '#!/bin/bash 
# 
#SBATCH --job-name=T-'$T'
#SBATCH --output=out 
#SBATCH -e errorLog 
#SBATCH -n 16 
#SBATCH -p short
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
#SBATCH -t 1:00:00 

module load OpenMPI 
mpirun ./h2dmmulti 3200 > log' > $DIRNAME/T-$T.job

done

