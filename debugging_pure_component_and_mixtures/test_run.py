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
import sys
import argparse

# No, I know this isn't the best way to do this, but between having my runnign script in a separate directory along with not quite having all my needed variables set using import, this seems like the best way to handle it right now. I blame Python for not working like other high level scripting languages
exec(open("gcmc.py").read())


# Default Thermo Settings
yco = 0.5
P_res = 200*10**5 #[Pa]
T = 45 + 273.15 #K

# Simulation Settings
mega_verbose = False
verbose = False
s_box = 57.15
N_max = 50000
N_moves = 1
N_equil = 10
N_prod = 10 #int( np.round( 500000/N_moves) )

# Carbon Wall Settings
sf = False
del_sf = 3.35 #[A]
rho_sf = 0.114 #[A^-3]
W = 5*3.8 #[A] relative to diameter of methane 3.80A


parser = argparse.ArgumentParser()
parser.add_argument("-y", "--yco", action = "store", dest="yco",type=float)
parser.add_argument("-p", "--pres", action = "store", dest="P_res",type=float)
parser.add_argument("-t", "--tres", action = "store", dest="T",type=float)
parser.add_argument("-v", "--verbose", action = "store_true",dest="v",default=False)
parser.add_argument("-m", "--mverbose", action = "store_true",dest="m",default=False)
args = parser.parse_args()

if args.yco:
    yco = args.yco
if args.P_res:
    P_res = args.P_res*10**5
if args.T:
    T = args.T + 273.15
if args.v:
    verbose = True
if args.m:
    mega_verbose = True

print( "yco", yco )
print( "P_res", P_res )
print( "T", T )

# Derived Values & Constants
Vol = s_box**3
fco, fme = PR_Fugacity( P_res/10**5, T, yco )
fco = fco*10**5
fme = fme*10**5
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


rhocov,rhomev,Env,Pv,Ncov, Nmev = mc_run(verbose = verbose)


