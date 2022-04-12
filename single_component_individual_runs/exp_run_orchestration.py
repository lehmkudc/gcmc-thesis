import subprocess
import pandas as pd
#import numpy as np
#import sys
#import os
from multiprocessing import Process


exp_filepath = "exp_list.csv"

def run_single( exp_list, i):
    yco = exp_list.yco[i]
    p_mpa = exp_list.p_mpa[i] 
    t_c = exp_list.t_c[i]
    s_box = exp_list.s_box[i]
    n_moves = exp_list.n_moves[i]
    n_equil = exp_list.n_equil[i]
    n_prod = exp_list.n_prod[i]
    rep = exp_list.rep[i]
    filepath = "data/" + str(exp_list.exp[i]) + ".csv"
    
    shellString = ("python ../one_run.py -y " + str(yco) +
                   " -p " + str(p_mpa) + " -t " + str(t_c) + 
                   " -s " + str(s_box) + " -m " + str(n_moves) + 
                   " -e " + str(n_equil) + " -o " + str(n_prod) + 
                   " -r " + str(rep) + " -f " + str(filepath)
                   )
    print( shellString )
    output = subprocess.check_output(shellString, shell = True)
    exp_list.at[i,"complete"] = 1
    print( output )
    
    return( exp_list )
    #return( output )

# if __name__ == '__main__':
#     experiments = pd.read_csv( "exp_list.csv")
    
#     process_list = []
#     for i in range(experiments.shape[0]):
#         p =  Process(target= run_single, args = (experiments,i))
#         p.start()
#         process_list.append(p)
    
#     for process in process_list:
#         process.join()
 
experiments = pd.read_csv(exp_filepath)
       
for i in range(experiments.shape[0]):
    
    if experiments.complete[i] == 0:
        experiments = run_single( experiments, i)
        experiments.to_csv(exp_filepath, index = False)
        
        
print( experiments )