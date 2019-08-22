#!/bin/bash

### this script reads the binary measurement files in a given T range 
### and produces plain text files with all data needed for python plotting

if [ "$#" -ne 3 ]; then
    echo "Usage: <script name> Tmin Tmax dT"
	exit 1
fi

TEMPDIR='temp'

if [ ! -d "$TEMPDIR" ]; then
  mkdir $TEMPDIR
fi

LABELS_H='T ?? Avg(hsq) error autoc ind.meas betaG2 g3 RGscale3d gsq g1sq ytsq Lb Lf'


Tmin=$1
Tmax=$2
dT=$3

for T in $(seq $Tmin $dT $Tmax);
do

	MEASUREMENTDIR='T-'$T

	if [ ! -d "$MEASUREMENTDIR" ]; then
  		continue
	fi

	aa -d4 $MEASUREMENTDIR/measure > $TEMPDIR/h1temp_$T.dat
	aa -d7 $MEASUREMENTDIR/measure > $TEMPDIR/h2temp_$T.dat
	
	awk '{print}' $TEMPDIR/h1temp_$T.dat >> $TEMPDIR/data_h1.dat
	awk '{print}' $TEMPDIR/h2temp_$T.dat >> $TEMPDIR/data_h2.dat

	echo $T >> $TEMPDIR/T.dat

	# get betaG2, g3 for removing lattice divergence,
	# as well as 4d gsq, g1sq, ytsq, Lb, Lf for relating to 4d fields 
	awk 'NR==1 {print $2}' beta >> $TEMPDIR/tempBeta
	awk -v T="$T" '$1 + 0 == T { print $2 "\t" $11}' input_couplings >> $TEMPDIR/input_temp0

done

paste $TEMPDIR/tempBeta $TEMPDIR/input_temp0 > $TEMPDIR/input_temp
rm $TEMPDIR/tempBeta
rm $TEMPDIR/input_temp0

paste $TEMPDIR/T.dat $TEMPDIR/data_h1.dat > h1_T.dat
paste $TEMPDIR/T.dat $TEMPDIR/data_h2.dat > h2_T.dat

echo $LABELS_H > labels_h

# delete unnecessary temp files 
rm $TEMPDIR/h1temp*.dat
rm $TEMPDIR/h2temp*.dat
rm $TEMPDIR/data_h1.dat
rm $TEMPDIR/data_h2.dat
rm $TEMPDIR/T.dat

# add data for input betaG, g3, RG scale to the final data files.
# (they are needed for removing lattice divergence)
paste h1_T.dat $TEMPDIR/input_temp > h1_T_temp.dat
paste h2_T.dat $TEMPDIR/input_temp > h2_T_temp.dat
mv h1_T_temp.dat h1_T.dat
mv h2_T_temp.dat h2_T.dat

# problem: sometimes the autocorrelation column comes with numbers in scientific notation:
# this messes up the column structure. convert to decimal:
sed -i 's/(/( /g' h1_T.dat 
sed -i 's/(/( /g' h2_T.dat 
awk '{$6=sprintf("%.2f", $6)")"}1' h1_T.dat > h1_Temp.dat
awk '{$6=sprintf("%.2f", $6)")"}1' h2_T.dat > h2_Temp.dat

mv h1_Temp.dat h1_T.dat
mv h2_Temp.dat h2_T.dat

# run the python plotting script
#python plot_continuum.py














