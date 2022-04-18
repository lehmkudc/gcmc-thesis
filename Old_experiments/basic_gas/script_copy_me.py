from random import random, seed, randrange
from math import floor, pi
import matplotlib.pyplot as plt
import numpy as np
import os
from time import time
import pandas as pd
from scipy.integrate import simps, trapz, cumtrapz
from scipy.ndimage.filters import gaussian_filter1d
from time import time

string = os.path.basename(__file__)
string = string.replace(".py", "")
print( string )

exec(open("../carbon_walls/gcmc.py").read())

# filename = "../data/" + "" + ".csv"
data_dump = pd.read_csv(filename).copy()

for i in range(data_dump.shape[0]):
    if np.isnan( data_dump.yco_sim[i]):
        yco = data_dump.yco_res[i]
        P_res = data_dump.P_res[i]*10**5 #[Pa]
        T = data_dump.T_res[i] + 273.15 #K
        fco, fme = PR_Fugacity( P_res/10**5, T, yco )
        fco = fco*10**5
        fme = fme*10**5
        
        sf = False
        mega_verbose = False

        s_box = 57.15
        N_max = 50000
        Vol = s_box**3
        kb = 1.3806*10**(7) #[Pa*A^3/K]
        Nco = 0 #floor(fco*Vol/kb/T)
        Nme = 0 #floor(fme*Vol/kb/T)
        Nc = 0
        rc = s_box
        beta = 1/T
        zz_co = beta*fco
        zz_me = beta*fme
        delta = 1
        pi_move = 0.5

        N_moves = 100
        if T == 30+273.15:
            N_equil = int( np.round( 400000/N_moves) )
        if T == 45+273.15:
            N_equil = int( np.round( 350000/N_moves) )
        if T == 60+273.15:
            N_equil = int( np.round( 300000/N_moves) )
        N_prod = int( np.round( 10000/N_moves) )

        t0 = time()
        rhocov,rhomev,Env,Pv,Ncov, Nmev = mc_run()
        t = time()
        
        Nco = Ncov.mean()
        Nme = Nmev.mean()
        data_dump.yco_sim[i] = Nco/(Nco + Nme)
        data_dump.P_sim[i] = Pv.mean()*10 #[bar]
        data_dump.E_sim[i] = Env.mean()
        data_dump.rhoco[i] = rhocov.mean()
        data_dump.rhome[i] = rhomev.mean()
        data_dump.yco_s[i] = (rhocov/(rhocov + rhomev)).std()
        data_dump.P_s[i] = Pv.std()*10 #[bar]
        data_dump.rhoco_s[i] = rhocov.std()
        data_dump.rhome_s[i] = rhomev.std()
        data_dump.E_s[i] = Env.std()
        data_dump.time[i] = t - t0
        print( i, data_dump.yco_sim[i],  data_dump.P_sim[i], data_dump.rhoco[i], data_dump.rhome[i], data_dump.time[i] )
        #filename = "../data_output/yco" + str(yco) + "P" + str(P) + "/run" + str(i) + ".csv"
        data_dump.to_csv( filename, index = False)