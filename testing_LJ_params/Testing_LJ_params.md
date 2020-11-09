# Investigating LennardJones Parameters for CO2

During the course of my previous experiments, I've noticed a serious discrepancy between expected data and simulated results for CO2. After testing another LJ set of values, the model behaved better but not particularly well.
I'm running some comparisons to see
a) What the best LJ parameters are for CO2 that I ahve access to
b) How sensitive to LJ parameters this simulation is.

Ideally a pure component CO2 should behave very similarly to NIST thermophysical data.


Steps of Experiment: 
1. Determine the number of iterations needed to equilibrate pure component CO2 at max pressure
2. Run plenty of simulations at a range of pressures in triplicate
3. Compare the resulting isotherms to eachother as well as NIST Standard Data.