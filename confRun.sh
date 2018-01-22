#!/bin/bash
echo ""
echo Set up model parameters and calculate terms of the series expoantion 
echo ""

# Compile the source file
g++ -std=c++11 -o NumericalAnalysis/calculateTerms.out NumericalAnalysis/calculateTerms.cpp -O2 || exit 1

# muMax sets the limit of the service rate range
muMax=10.0; echo Range of service rates is [0...$muMax]

# carRate sets average rate of new car arrival
carRate=0.1; echo average rate of new car arrival is $carRate

# lam defiles arrival rate
lam=3.0

# K defines the maximum number of cars in the system
K=3; echo Max number of cars in the system is $K

# C defines capacity of a single queue
C=3; echo Queueing capacity is $C

# M defines number of zones
M=3; echo Number of zones is $M

# nTerms defines number of terms to calculate
nTerms=30; echo $nTerms terms will be calculated

# stopCrit defines stopping criteria
stopCrit=0.00001; echo Stopping criteria is $stopCrit #0.000001

# w defines parameter of weighted jacoby
w=0.9; echo method: weighted Jacoby with parameter $w

### Simulaiton parameters
nSimSamples=20; echo nSimSamples = $nSimSamples defines number of samples
nSimLoops=3000000; echo nSimLoops = $nSimLoops defines number of simulation loops to obtain one sample

# outFile is a file where output is written. If empty, generated automaticly in the script
# outFile=testResult.dat; echo Output is written to $outFile in *.mat formal.
# To read results start matlab and load(outFile)
outFile=""
outFileSim=""

## Scenarios
# Optimistic DSRC 1/10/1 Mbits/sec, 400/1200/400 m, speed=90 km/h
# results in beta = [1/16, 1/48, 1/16]
beta1=0.0625 ; beta2=0.0208; beta3=0.0625
echo Rates of switching zones beta = [$beta1, $beta2, $beta3]
mbits1=1
mbits2=10
mbits3=1
echo Transmission rate in zones [$mbits1, $mbits2, $mbits3]

# Pessimistic DSRC. 1/3/1 mbits, 150/300/150 m.

### Run numerical analysis
#/NumericalAnalysis/calculateTerms.out $muMax $carRate\
#				$K $C $M $lam $nTerms  \
#				$stopCrit $w  \
#				$beta1 $beta2 $beta3 \
#				$mbits1 $mbits2 $mbits3\
#				$outFile

# Run simulation
nohup matlab -nodisplay -r \ "addpath('Simulation'); run_simulation(\
				'$muMax', '$carRate',\
				'$K', '$C', '$M', '$lam', '$beta1', '$beta2', '$beta3', \
				'$mbits1', '$mbits2', '$mbits3', \
				'$nSimSamples', '$nSimLoops', '$outFileSim')"	> ./simlog.txt </dev/null &				
				
#nohup ./calculateTerms.out $muMax $K $C $M $lam $nTerms $outFile $stopCrit $w $beta $>& tail.txt &






