#!/bin/bash

# this script prepares a weightfile for kari's multicanonical simulation

if [ "$#" -ne 3 ]; then
    echo "Usage: <script name> min max bins"
	exit 1
fi

min=$1
max=$2
bins=$3

if [ -f weight ]; then
    rm weight
fi

if [ -f temp ]; then
    rm temp
fi

# calculate width of each bin.
width=`echo "scale=10; ($max - $min)/$bins" | bc`


count=-1 ## for my code, write number of bins in the range; for Kari write bins+1

for i in `seq $min $width $max`
do
	echo $i    0 >> temp
	((count++))
done


delta=0.5
last_max=0

sed "1 i$count $delta $last_max $min $max" temp > weight  

rm temp
