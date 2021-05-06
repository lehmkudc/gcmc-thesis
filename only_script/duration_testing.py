from random import random, seed, randrange
from math import floor, pi
import numpy as np
import os
from time import time
import pandas as pd
from time import time
import sys
import argparse

# No, I know this isn't the best way to do this, but between having my runnign script in a separate directory along with not quite having all my needed variables set using import, this seems like the best way to handle it right now. I blame Python for not working like other high level scripting languages
exec(open("gcmc.py").read())

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--rep", action = "store", dest="rep",type=float)
args = parser.parse_args()

rep = args.rep


yco = 1
P_res = 200*10**5 #[Pa]
T = 30 + 273.15 #K
fco, fme = PR_Fugacity( P_res/10**5, T, yco )
fco = fco*10**5
fme = fme*10**5

del_sf = 3.35 #[A]
rho_sf = 0.114 #[A^-3]
W = 5*3.8 #[A] relative to diameter of methane 3.80A
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

N_moves = 1000
N_equil = 0
N_prod = int( np.round( 2000000/N_moves) )


rhocov,rhomev,Env,Pv,Ncov, Nmev = mc_run(verbose = True)

output = pd.DataFrame()
output['rhocov'] = rhocov
output['rhomev'] = rhomev
output['Env'] = Env
output['Pv'] = Pv
output['Ncov'] = Ncov
output['Nmev'] = Nmev

output.to_csv("duration_test" + str(rep) + ".csv", index = False)