#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: <script name> <input file> <T> <beta_G>"
	exit 1
fi

inputfile=$1
T=$2
beta=$3

read_params.py $inputfile $T $beta

# read the file created by the python script

latticefile='params_lattice.dat'

params=$(tail -n +2 $latticefile)

T=$(echo $params | awk {'print $1'})
beta=$(echo $params | awk {'print $2'})
g=$(echo $params | awk {'print $3'})
gp=$(echo $params | awk {'print $4'})
muphisq=$(echo $params | awk {'print $5'})
muSigmasq=$(echo $params | awk {'print $6'})
lam=$(echo $params | awk {'print $7'})
b4=$(echo $params | awk {'print $8'})
a2=$(echo $params | awk {'print $9'})
a=$(echo $params | awk {'print $10'})


# write these into the config file (assuming it exists)
cp config config_bu

## this replaces the line after "gauge coupling" with g $g 
sed -i -e "/gauge couplings/!b;n;c\betasu2 $beta" config

sed -i -e "/# doublet/!b;n;c\msq $muphisq" config
sed -i -e "/# doublet/!b;n;n;c\lambda $lam" config

sed -i -e "/# triplet/!b;n;c\msq_triplet $muSigmasq" config
sed -i -e "/# triplet/!b;n;n;c\b4 $b4" config
sed -i -e "/# triplet/!b;n;n;n;c\a2 $a2" config

