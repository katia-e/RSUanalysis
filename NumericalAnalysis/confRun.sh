#!/bin/bash
echo ""
echo Set up model parameters and calculate terms of the series expoantion 
echo ""

# Compile the source file
g++ -std=c++11 -o calculateTerms.out calculateTerms.cpp -O3

# muMax sets the limit of the service rate range
muMax=10.0; echo Range of service rates is [0...$muMax]

# K defines the maximum number of cars in the system
K=10; echo Max number of cars in the system is $K

# C defines capacity of a single queue
C=2; echo Queueing capacity is $C

# M defines number of zones
M=3; echo Number of zones is $M

# outFile is a file where output is written
outFile=testTesult.dat; echo Output is written to $outFile

# nTerms defines number of terms to calculate
nTerms=30; echo $nTerms terms will be calculated

# stopCrit defines stopping criteria
stopCrit=0.000001; echo Stopping criteria is $stopCrit

# w defines parameter of weighted jacoby
w=0.9; echo method: weighted Jacoby with parameter $w

nohup ./calculateTerms.out $muMax $K $C $M 3.0 $outFile $nTerms $stopCrit $w >& tail.txt &






