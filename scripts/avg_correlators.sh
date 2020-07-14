#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: <script name> <input file>"
	exit 1
fi

fname=$1

## cut out identifier lines
grep -v 'Measurement id' $fname > temp 

## remove blank lines
sed -i '/^$/d' temp

## sort based on the value in first column (should be distance d)
sort -k1 -n temp > temp1 

## find largest distance, which is now at the end of temp1
max=$(tail -n1 temp1 | awk '{print $1}')

rm temp1

## write header for results file
echo 'd Sigma^2 Sigma_err photon_re re_err photon_im im_err' > correl_avg


## now for each d, calculate avg of correlators
for (( d=0; d<=$max; d++ ))
do

	## find all lines with distance=d 
	awk -v d="$d" '{if($1==d) print}' temp > corr_temp

	## averages using aa. columns 2 and 3 of the output are avg and error
	# Sigma correlator
	sig=$(aa -d2 -q corr_temp)
	sig=$(echo $sig | awk '{print $2, $3}')

	
	# photon correlator, Tr Sigma F_ij, real and imag parts
	re=$(aa -d3 -q corr_temp)
	re=$(echo $re | awk '{print $2, $3}')
	im=$(aa -d4 -q corr_temp)
	im=$(echo $im | awk '{print $2, $3}')
	

	# print to results file
	echo "$d $sig $re $im" >> correl_avg
	
	rm corr_temp

done 


rm temp

