# RSUanalysis
Numerical performance evaluation of Drive-Thru Internet scenario

Run analysis and simulation with bash script:
``` bash
./confRun.sh
```

All input parameters can be set inside `./confRun.sh`
### Changable parameters:
 * `car_rate` - rate of the vehicles entering the system
 * `K` - maximum number of castomers in the system
 * `C` - buffer capacity
 * `M` - number of zones
 * `lamJac` - terms are calculated in the region of this value
 * `nTerms` - number of terms to calculate
 * `stopCrit` - stopping criteria for terms calculation 
 * `w` - parameter of weighted jacoby
 
 Simulation:
 * `lamMax` - maximum Lambda for simulation
 * `nSimSamples` - number of samples in simulation
 * `nSimLoops` - number of simulation loops
