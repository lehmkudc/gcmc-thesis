# Experiment 3: Testing Efficacy of Simulation for Pure Component Gas

## Premise
In Experiment 2, I was uncomfortable with the significant bias shown by simulating both gas mixtures and single-component gases which I recall being much closer in previous iterations of my code. In part to interrogate any bugs, thermodynamic model mistakes, or improper model parameters, I'll be focusing on pure component gases compared to PR equation of state, as well as NIST Thermophysical property data.

An "Ideal" result would be no systematic bias between my simulation and the results given by NIST and Peng-Robinson. The variability between results would decrease given simulation complexity and computation time. This would prove that my simulation code can properly represetnt a non-adsorption system accurately. 

## Goals
1. Determine if my simulation code consistently deviates from accepted thermophysical properties for Pure CO2 and Pure Methane.
2. Characterize these differences using a series of simulation integrity tests.
3. Compare results across a series of simulation conditions and thermodynamic state properties.

## Sub-Experiments:
1. Visually determine appropriate number of equilibration steps through testing at Sim extremes
2. Give a baseline first-pass result with "fast" simulation results
3. Analyze compared to NIST & Theoretical Data.
4. Increase the Box Size
5. Increase the Equilibration Duration

## Conditions
* Mole Fractions: 0, 1
* Pressures: 5, 10, 15, ..., 200 MPa
* Temperatures: 30, 45, 60 C
* Steps Per Cycle: 1K
* Equilibraiton Steps: 1M, 1.5M, 2M (Determined by Sub Experiment)
* Production Cycles: 100, 1K, 10K
* Unit Cell Length: 40A, 45A, 50A

