#!/bin/bash

### this script reads the binary measurement files in a given T range 
### and produces plain text files with all data needed for python plotting

if [ "$#" -ne 2 ]; then
    echo "Usage: <script name> Tmin Tmax"
	exit 1
fi

TEMPDIR='temp'

if [ ! -d "$TEMPDIR" ]; then
  mkdir $TEMPDIR
fi

LABELS_RES='T phisq phisq_err S S_err S2 S2_err beta RGscale gsq gpsq mphisq lam mSsq b1 b3 b4 a1 a2'

RESFILE='condensates_T.dat'
SUSCFILE='susc_T.dat'

echo $LABELS_RES > $RESFILE
echo $LABELS_RES > $SUSCFILE

Tmin=$1
Tmax=$2

## Data columns
col_phisq=7
col_S=9
col_Ssq=10


for dir in */
do

	dir=${dir%*/} # remove trailing '/'

	measfile=$dir/measure
	paramfile=$dir/params_continuum.dat

	if [[ ! -f "$measfile" || ! -f "$paramfile" ]]; then
		echo "Skipping directory $dir"
  		continue
	fi
	
	## read continuum parameters from params_continuum.dat ; need these later to subtract divergence
	T=$(grep T $paramfile | awk '{print $2}')
	RGscale=$(grep RGscale $paramfile | awk '{print $2}')
	gsq=$(grep gsq $paramfile | awk '{print $2}')
	gpsq=$(grep gpsq $paramfile | awk '{print $2}')
	mphisq=$(grep mphisq $paramfile | awk '{print $2}')
	mSsq=$(grep mSsq $paramfile | awk '{print $2}')
	lambda=$(grep lambda $paramfile | awk '{print $2}')
	b1=$(grep b1 $paramfile | awk '{print $2}')
	b3=$(grep b3 $paramfile | awk '{print $2}')
	b4=$(grep b4 $paramfile | awk '{print $2}')
	a1=$(grep a1 $paramfile | awk '{print $2}')
	a2=$(grep a2 $paramfile | awk '{print $2}')

	## also beta
	beta=$(grep betasu2 $dir/config | awk '{print $2}')
	
	
	if (( $(echo "$T > $Tmax" |bc -l) || $(echo "$T < $Tmin" |bc -l) )); then
		continue
	fi


	### get volume averages, their susceptibilities and errors

	## <phisq>
	phisq=$(aa -d$col_phisq $measfile | awk '{print $2,$3}')
	phisq_susc=$(aa -d$col_phisq -2 $measfile | awk '{print $2,$3}')
	## <S> 
	S=$(aa -d$col_S $measfile | awk '{print $2,$3}')
	S_susc=$(aa -d$col_S -2 $measfile | awk '{print $2,$3}')
	## <S^2>
	S2=$(aa -d$col_Ssq $measfile | awk '{print $2,$3}')
	S2_susc=$(aa -d$col_Ssq -2 $measfile | awk '{print $2,$3}')


	## finally, append to the results files
	echo "$T $phisq $S $S2 $beta $RGscale $gsq $gpsq $mphisq $lambda $mSsq $b1 $b3 $b4 $a1 $a2" >> $RESFILE
	echo "$T $phisq_susc $S_susc $S2_susc $beta $RGscale $gsq $gpsq $mphisq $lambda $mSsq $b1 $b3 $b4 $a1 $a2" >> $SUSCFILE

done










