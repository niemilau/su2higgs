#!/bin/bash

# this script prepares a weightfile for multicanonical simulation

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


count=0

for i in `seq $min $width $max`
do
	echo $i    0 >> temp
	((count++))
done

sed "1 i$count" temp > weight  

rm temp

