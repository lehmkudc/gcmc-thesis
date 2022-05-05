import argparse
import pandas as pd
import numpy as np
exec(open("gcmc.py").read())

parser = argparse.ArgumentParser()
parser.add_argument("-y", "--yco", action = "store", dest="yco",type=float, default = 1)
parser.add_argument("-p", "--p_bar", action = "store", dest="p_bar",type=float, default = 20)
parser.add_argument("-t", "--t_c", action = "store", dest="t_c",type=float, default = 45)
parser.add_argument("-s", "--s_box", action = "store", dest="s_box",type=float, default = 34)
parser.add_argument("-m", "--n_moves", action = "store", dest="n_moves",type=int, default = 1000)
parser.add_argument("-e", "--n_equil", action = "store", dest="n_equil",type=int, default = 0)
parser.add_argument("-o", "--n_prod", action = "store", dest="n_prod",type=int, default = 100000)
parser.add_argument("-r", "--rep", action = "store", dest="rep",type=int, default = 0)
parser.add_argument("-f", "--filepath", action = "store", dest="filepath",type=str, default="default.csv")
args = parser.parse_args()

yco = args.yco
P_bar = args.p_bar
T_C = args.t_c
s_box = args.s_box
n_moves = args.n_moves
n_equil = args.n_equil
n_prod = args.n_prod
rep = args.rep
filepath = args.filepath



if filepath == "default":
    filepath = ("data/" + "yc0" + str(yco) + "_pbar" + str(P_bar) + 
                "_tc" + str(T_C) + "_sbox" + str(s_box) + "_nm" + str(n_moves) + 
                "_ne" + str(n_equil) + "_np" + str(n_prod) + "_r" + str(rep) +
                ".csv")
    
    
    
print( filepath )

P_res = P_bar*10**5 #[Pa]
T = T_C + 273.15 #K
fco, fme = PR_Fugacity( P_res/10**5, T, yco )
fco = fco*10**5
fme = fme*10**5

del_sf = "Nothing"#3.35 #[A]
rho_sf = "Nothing"#0.114 #[A^-3]
W = "Nothing"#5*3.8 #[A] relative to diameter of methane 3.80A
sf = False
mega_verbose = False

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


rhocov,rhomev,Env,Pv,Ncov, Nmev = mc_run(verbose = True)

output = pd.DataFrame()
output['rhocov'] = rhocov
output['rhomev'] = rhomev
output['Env'] = Env
output['Pv'] = Pv
output['Ncov'] = Ncov
output['Nmev'] = Nmev

print( output )

output.to_csv( filepath, index = False)