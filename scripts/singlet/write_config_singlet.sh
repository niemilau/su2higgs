#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: <script name> <input file> <T> <beta_G>"
	exit 1
fi

inputfile=$1
T=$2
beta=$3

./params_singlet.py $inputfile $T $beta

# read the file created by the python script

latticefile='params_lattice.dat'

T=$(grep T $latticefile | awk '{print $2}')
a=$(grep spacing $latticefile | awk '{print $2}')
gsq=$(grep gsq $latticefile | awk '{print $2}')
gpsq=$(grep gpsq $latticefile | awk '{print $2}')
mphisq=$(grep mphisq $latticefile | awk '{print $2}')
mSsq=$(grep mSsq $latticefile | awk '{print $2}')
lambda=$(grep lambda $latticefile | awk '{print $2}')
b1=$(grep b1 $latticefile | awk '{print $2}')
b3=$(grep b3 $latticefile | awk '{print $2}')
b4=$(grep b4 $latticefile | awk '{print $2}')
a1=$(grep a1 $latticefile | awk '{print $2}')
a2=$(grep a2 $latticefile | awk '{print $2}')


# write these into the config file (assuming it exists)
cp config config_bu

sed -i "s/at T = .*/at T = $T, a = $a/g" config
sed -i "s/msq .*/msq $mphisq/g" config
sed -i "s/lambda .*/lambda $lambda/g" config
sed -i "s/msq_s .*/msq_s $mSsq/g" config
sed -i "s/b1_s .*\+/b1_s $b1/g" config
sed -i "s/b3_s .*\+/b3_s $b3/g" config
sed -i "s/b4_s .*\+/b4_s $b4/g" config
sed -i "s/a1_s .*\+/a1_s $a1/g" config
sed -i "s/a2_s .*\+/a2_s $a2/g" config

