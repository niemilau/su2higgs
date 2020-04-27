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
## my code: last bin ENDS at $max, so its start value is at $max-$width
lastbin=`echo "scale=10; ($max-$width)" | bc`
## while for Kari's code the last bin starts at max

for i in `seq $min $width $lastbin`
do
	echo $i    0 >> temp
	((count++))
done

# on the first line write number of bins, increment factor, last_max, min, max, min_abs, max_abs
sed "1 i$count 0.5 0 $min $max $min $max" temp > weight  
## for Kari's code write instead number of bins+1 and nothing else

rm temp

