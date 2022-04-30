import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

dataPath = "single_component_individual_runs/longer_tests/data/"

allData = pd.DataFrame()

for dataFile in os.listdir(dataPath ):
    df = pd.read_csv( os.path.join( dataPath, dataFile ) )
    df['source'] = dataFile
    df['exp'] = int( dataFile.replace( ".csv", "" ) )
    df['step'] = range( df.shape[0] )
    
    allData = allData.append( df )
    
avgData = allData.groupby(['step']).agg({ 'rhocov': ['mean', 'std'], 'Env': ['mean', 'std'], 'Pv': ['mean', 'std']})
avgData.columns = ['rhocov_mean','rhocov_std','Env_mean','Env_std','Pv_mean','Pv_std']
avgData = avgData.reset_index()


nist_co2 = pd.read_csv( 'reference_data/c02.txt', delimiter = '\t')

nist_co2.head()
nist_co2.columns

target_co2 = nist_co2[ nist_co2['Pressure (bar)'] == 5.0]

nist_rhoco = target_co2['Density (mol/l)']
nist_mol = nist_rhoco
nist_U = target_co2['Internal Energy (J/mol*K)']


#### Density ####
plt.figure()
for run in sorted( allData.exp.unique() ):
    df = allData[allData.exp == run ]
    plt.plot(
        df.step,
        df.rhocov,
        alpha = 0.6,
        linewidth = 1,
        color = "black"
    )
    
plt.plot( avgData.step, avgData.rhocov_mean, alpha = 0.9, linewidth = 2 )

plt.title("Average Density of Pure CO2")
plt.xlabel("Iterations [10^4]")
plt.ylabel("Density [Parts Per A^3]")
plt.show()


#### Internal Energy ####
plt.figure()
for run in sorted( allData.exp.unique() ):
    df = allData[allData.exp == run ]
    plt.plot(
        df.step,
        df.Env,
        alpha = 0.6,
        linewidth = 1
    )


plt.plot(
    avgData.step, avgData.Env_mean, alpha = 0.9, linewidth = 2
)

plt.title("Average Internal Energy")
plt.xlabel("Iterations [10^4]")
plt.ylabel("Internal Energy [J/K]")
plt.show()



#### Internal Energy ####
plt.figure()
for run in sorted( allData.exp.unique() ):
    df = allData[allData.exp == run ]
    plt.plot(
        df.step,
        df.Pv,
        alpha = 0.6,
        linewidth = 1,color="black"
    )


plt.plot(avgData.step, avgData.Pv_mean, alpha = 0.9, linewidth = 2, color = "limegreen")
plt.plot( avgData.step, np.full(100,5), alpha = 0.9, linewidth = 2, color = "red" )

plt.title("Average Pressure")
plt.xlabel("Iterations [10^4]")
plt.ylabel("Pressure [MPa]")
plt.show()



