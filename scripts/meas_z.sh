#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: <script name> <input file> <measurement id>"
	exit 1
fi

fname=$1
id=$2

if [ ! -f "$fname" ]; then
		echo "No such file!!"
		exit 2
fi

## full header for the measurement 
header=$(echo "Measurement $2:")

## line number of the matching header? (starts from 1)
start_line=$(grep -n "$header" $fname | cut -d: -f1)


## if there were more than one match, return error. Mixed measurements from different simulations?
[[ $start_line == *$'\n'* ]] && echo "Error! Found more multiple lines matching header \"$header\"" && exit 3

[[ $start_line == "" ]] && echo "No such measurement!" && exit 4

start_line=$(($start_line+1)) ## don't print the header

## next measurement?
end=$(($id + 1))
header=$(echo "Measurement $end:")

end_line=$(grep -n "$header" $fname | cut -d: -f1)

## if the next measurement is not found, we prolly asked for the last measurement. just print till EOF
[[ $end_line == "" ]] && tail -n +$start_line $fname


## else, print lines between start and end lines
awk -v s="$start_line" -v e="$end_line" 'NR>=s&&NR<e' $fname




