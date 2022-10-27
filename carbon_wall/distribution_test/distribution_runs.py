import argparse
import pandas as pd
import numpy as np
from time import time
exec(open("gcmc.py").read())


parser = argparse.ArgumentParser()
parser.add_argument("-y", "--yco", action = "store", dest="yco",type=float, default = 0.5)
parser.add_argument("-p", "--p_bar", action = "store", dest="p_bar",type=float, default = 20)
parser.add_argument("-t", "--t_c", action = "store", dest="t_c",type=float, default = 45)
parser.add_argument("-s", "--s_box", action = "store", dest="s_box",type=float, default = 50)
parser.add_argument("-m", "--n_moves", action = "store", dest="n_moves",type=int, default = 10000)
parser.add_argument("-e", "--n_equil", action = "store", dest="n_equil",type=int, default = 5000000)
parser.add_argument("-o", "--n_prod", action = "store", dest="n_prod",type=int, default = 1000000)
parser.add_argument("-r", "--row", action = "store", dest="row",type=int, default = 0)
parser.add_argument("-f", "--filepath", action = "store", dest="filepath",type=str, default="default")
parser.add_argument("-c", "--sf", action = "store", dest="sf",type=bool, default = False)
parser.add_argument("-d", "--del_sf", action = "store", dest="del_sf",type=float, default = 3.35)
parser.add_argument("--rho_sf", action = "store", dest="rho_sf",type=float, default = 0.114)
parser.add_argument("-w", "--w", action = "store", dest="w",type=float, default = 19)
args = parser.parse_args()


yco = args.yco
P_bar = args.p_bar
T_C = args.t_c
s_box = args.s_box
n_moves = args.n_moves
n_equil = args.n_equil
n_prod = args.n_prod
row = args.row
filepath = args.filepath

P_res = P_bar*10**5 #[Pa]
T = T_C + 273.15 #K
fco, fme = PR_Fugacity( P_res/10**5, T, yco )
fco = fco*10**5
fme = fme*10**5

mega_verbose = False

sf = args.sf
del_sf = args.del_sf #[A]
rho_sf = args.rho_sf #[A^-3]
W = args.w #[A] relative to diameter of methane 3.80A


e_c = 28.0 # eps over kb[K] # Activated is 89.44
s_c = 3.4  # sigma [A]

# Surface of carbons
e_csf = 28.0 # es over kb[K]
s_csf = 3.40 # sigma [A]

e_meco = np.sqrt(e_me*e_co)
s_meco = (s_me+s_co)/2

e_cme = np.sqrt(e_me*e_c)
s_cme = (s_me+s_c)/2

e_cco = np.sqrt(e_co*e_c)
s_cco = (s_co + s_c)/2

e_cosf = np.sqrt(e_co*e_csf)
s_cosf = (s_co + s_csf)/2

e_mesf = np.sqrt(e_me*e_csf)
s_mesf = (s_me + s_csf)/2

s_box = s_box #[A]
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

N_moves = n_moves
N_equil = int( np.round( n_equil/N_moves) )
N_prod = int( np.round( n_prod/N_moves) )

t0 = time()
rhocov,rhomev,Env,Pv,Ncov, Nmev = mc_run(verbose = True)

print( rhocov[0] )
tf = time()


output_co = pd.DataFrame()
output_co['Xco'] = np.array(Xco)[:Nco]
output_co['Yco'] = np.array(Yco)[:Nco]
output_co['Zco'] = np.array(Zco)[:Nco]

output_co.to_csv( filepath + "_co.csv", index = False )

output_me = pd.DataFrame()
output_me['Xme'] = np.array(Xme)[:Nme]
output_me['Yme'] = np.array(Yme)[:Nme]
output_me['Zme'] = np.array(Zme)[:Nme]

output_me.to_csv( filepath + "_co.csv", index = False )