import pandas as pd
import numpy as np
exec(open("../gcmc.py").read())

yco = [0,1]
t_c = [30,45,60]
p_bar = np.linspace(5,200,40)

i = 0
totalLength = len(yco)*len(t_c)*len(p_bar)
yv = [0]*totalLength
tv = [0]*totalLength
pv = [0]*totalLength 
zv = [0]*totalLength
rhov = [0]*totalLength

R = 0.0831446261815324 #[L bar/K mol]

for y in yco:
    for t in t_c:
        for p in p_bar:
            
            z = PR_Zmix( P_res= p, T = t+273.15, yco = y )
            yv[i] = y
            tv[i] = t
            pv[i] = p
            zv[i] = z
            rhov[i] = p/z/R/(t+273.15) #[mol/L]
            i = i+1
            
            
output = pd.DataFrame()
output['yco'] = yv
output['t_c'] = tv
output['p_bar'] = pv
output['z_mix'] = zv
output['rho'] = rhov

print( output )

output.to_csv("../reference_data/pr.csv")