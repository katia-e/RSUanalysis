#!/bin/bash
echo ""
echo Set up model parameters and calculate terms of the series expoantion 
echo ""

# Compile the source file
g++ -std=c++11 -o calculateTerms.out calculateTerms.cpp -O2

# muMax sets the limit of the service rate range
muMax=10.0; echo Range of service rates is [0...$muMax]

# lam defiles arrival rate
lam=3.0

# K defines the maximum number of cars in the system
K=10; echo Max number of cars in the system is $K

# C defines capacity of a single queue
C=3; echo Queueing capacity is $C

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

### Scenarios
# scenario 1. OptimisticDSRC 400/1200/400 m, speed=90 km/h
mbits=(1, 10, 1); echo Transmission rate in zones ${mbits[@]}
beta=(1/16, 1/48, 1/16); echo Rate of switching zones beta = ${beta[@]}

beta_path=scenarios/OptimisticDSRC.txt

# scenario 2. 1/3/1 mbits, 150/300/150 m.
# beta_path=PessimisticDSRC


./calculateTerms.out $muMax $K $C $M $lam $nTerms $outFile $stopCrit $w $beta
#nohup ./calculateTerms.out $muMax $K $C $M $lam $nTerms $outFile $stopCrit $w $beta $>& tail.txt &






