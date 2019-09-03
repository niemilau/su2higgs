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

LABELS_H='T avg error autoc ind.meas betaG gsq gpsq muphisq muSigmasq lambda b4 a2 RGscale'


Tmin=$1
Tmax=$2
dT=$3

# only searches directories in the exact T range; todo: improve this 
for T in $(seq $Tmin $dT $Tmax);
do

	MEASUREMENTDIR='T-'$T

	if [ ! -d "$MEASUREMENTDIR" ]; then
  		continue
	fi

	aa -d5 -s100 $MEASUREMENTDIR/measure > $TEMPDIR/h1temp_$T.dat
	aa -d8 -s100 $MEASUREMENTDIR/measure > $TEMPDIR/h2temp_$T.dat
	
	# cut aa label and write into a combined file
	awk '{print}' $TEMPDIR/h1temp_$T.dat | cut -d '1' -f2- >> $TEMPDIR/data_h1.dat
	awk '{print}' $TEMPDIR/h2temp_$T.dat | cut -d '1' -f2- >> $TEMPDIR/data_h2.dat

	echo $T >> $TEMPDIR/T.dat

	# get betaG and 3d MSbar parameters for subtracting lattice divergence
	awk 'NR>1 {print $2}' $MEASUREMENTDIR/params_lattice.dat >> $TEMPDIR/beta_temp
	awk 'NR>1 {print}' $MEASUREMENTDIR/params_MSbar.dat | cut -d ' ' -f2- >> $TEMPDIR/params_temp # cut T column

done

echo $LABELS_H > h1_T.dat 
echo $LABELS_H > h2_T.dat 

paste $TEMPDIR/T.dat $TEMPDIR/data_h1.dat $TEMPDIR/beta_temp $TEMPDIR/params_temp >> h1_T.dat
paste $TEMPDIR/T.dat $TEMPDIR/data_h2.dat $TEMPDIR/beta_temp $TEMPDIR/params_temp >> h2_T.dat

# delete unnecessary temp files 
rm $TEMPDIR/h1temp*.dat
rm $TEMPDIR/h2temp*.dat
rm $TEMPDIR/data_h1.dat
rm $TEMPDIR/data_h2.dat
rm $TEMPDIR/T.dat
rm $TEMPDIR/params_temp 
rm $TEMPDIR/beta_temp


# problem: sometimes the autocorrelation column comes with numbers in scientific notation:
# this messes up the column structure. convert to decimal:
#sed -i 's/(   /(/g' h1_T.dat   ## remove spaces after brackets
#sed -i 's/(   /(/g' h2_T.dat 
#awk 'NR>1 {$5=sprintf("%.2f", $5)")"}1' h1_T.dat > h1_Temp.dat
#awk 'NR>1 {$5=sprintf("%.2f", $5)")"}1' h2_T.dat > h2_Temp.dat

# just remove the bracket expression altogether
sed -i 's/(.*)/ /g' h1_T.dat
sed -i 's/(.*)/ /g' h2_T.dat










