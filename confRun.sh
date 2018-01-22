#!/bin/bash
echo ""
echo Set up model parameters and calculate terms of the series expoantion 
echo ""

# Compile the source file
g++ -std=c++11 -o NumericalAnalysis/calculateTermsLam.out NumericalAnalysis/calculateTermsLam.cpp -O2 || exit 1

# muMax sets the limit of the service rate range
lamJac=1.5; echo reference point for iterative method $lamJac

# carRate sets average rate of new car arrival
carRate=0.1; echo average rate of new car arrival is $carRate

# K defines the maximum number of cars in the system
K=5; echo Max number of cars in the system is $K

# C defines capacity of a single queue
C=3; echo Queueing capacity is $C

# M defines number of zones
M=3; echo Number of zones is $M

# nTerms defines number of terms to calculate
nTerms=10; echo $nTerms terms will be calculated

# stopCrit defines stopping criteria
stopCrit=0.0001; echo Stopping criteria is $stopCrit #0.000001

# w defines parameter of weighted jacoby
w=0.9; echo method: weighted Jacoby with parameter $w

### Simulaiton parameters
lamMax=10; echo Range of system loads [0...$lamMax]
nSimSamples=20; echo nSimSamples = $nSimSamples defines number of samples
nSimLoops=3000000; echo nSimLoops = $nSimLoops defines number of simulation loops to obtain one sample

## Scenarios
# Optimistic DSRC 1/10/1 Mbits/sec, 400/1200/400 m, speed=90 km/h => beta = [1/16, 1/48, 1/16]
beta1=0.0625 ; beta2=0.0208; beta3=0.0625
echo Rates of switching zones beta = [$beta1, $beta2, $beta3]
mbits1=1
mbits2=10
mbits3=1
echo Transmission rate in zones [$mbits1, $mbits2, $mbits3]

# Pessimistic DSRC. 1/3/1 mbits, 150/300/150 m.

### Run numerical analysis
./NumericalAnalysis/calculateTermsLam.out $lamJac $carRate\
				$K $C $M $lam $nTerms  \
				$stopCrit $w  \
				$beta1 $beta2 $beta3 \
				$mbits1 $mbits2 $mbits3

# Run simulation
#nohup matlab -nodisplay -r \ "addpath('Simulation'); run_simulation(\
#				'$lamMax', '$carRate',\
#				'$K', '$C', '$M', '$beta1', '$beta2', '$beta3', \
#				'$mbits1', '$mbits2', '$mbits3', \
#				'$nSimSamples', '$nSimLoops')"	> ./simlog.txt </dev/null &				
#				
#nohup ./calculateTerms.out $muMax $K $C $M $lam $nTerms $outFile $stopCrit $w $beta $>& tail.txt &






