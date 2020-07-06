# Experiment 1.5: Debugging a CO2 specific Simulation Discrepancy

For some reason when I ran some experiments my pure CO2 results were almost 2x the pressure than they were supposed to be, with densities to match. I'm testing some stuff using just pure components to see if there was just something fishy in my python environment or if there is a definite issue with my simulation.

Leads:
* This issue is not propogated with pure methane, so whatever is wrong is CO2 specific.
* The simulation code *was* mostly well behaved before, but there's a possibility that this behavior was masked because I wasn't running the simulation enough.