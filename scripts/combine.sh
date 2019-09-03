#!/bin/bash

## this script combines h1_T data files from different T ranges 

cat h1_*.dat > temph1.dat 
cat temph1.dat | sort > h1_T.dat
rm temph1.dat

cat h2_*.dat > temph2.dat 
cat temph2.dat | sort > h2_T.dat
rm temph2.dat
