#!/bin/bash

## this script goes through all directories labeled "T-something" and modifies the 
## max iteration count in config file. can be performed also during an ongoing simulation

if [ "$#" -ne 1 ]; then
    echo "Usage: <script name> <max iterations>"
	exit 1
fi

iters=$1


for dir in */
do 

	dir=${dir%*/} # dir is now the folder name without / at the end
	if [[ $dir != *"T-"* ]]; then
		continue
	fi
	
	if [ ! -f "$dir/config" ]; then
		echo "WARNING: Could not find file $dir/config" 
  		continue
	fi
 
	sed -i -e "/iterations/!b;n;c\iterations $iters" $dir/config

	## also set reset to 0 so that the iteration counter is not reset when restarting the simulation
	sed -i -e "/start new simulation /!b;n;c\reset 0" $dir/config

done

echo "Changed max iterations to $iters"
