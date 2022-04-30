import subprocess
import pandas as pd
#import numpy as np
#import sys
#import os
from multiprocessing import Pool, cpu_count


exp_filepath = "aspen_comparison/aspen_comparison.csv"
data_filepath = "aspen_comparison/data/"
run_script = "only_result.py"

# print( "Available CPU:", cpu_count() )

def run_single( exp_list, i):
    yco = exp_list.yco[i]
    p_bar = exp_list.p_bar[i] 
    t_c = exp_list.t_c[i]
    s_box = exp_list.s_box[i]
    n_moves = int( exp_list.n_moves[i] )
    n_equil = int( exp_list.n_equil[i] )
    n_prod = int( exp_list.n_prod[i] )
    filepath = data_filepath + "/" + str(exp_list.exp[i]) + ".csv"
    
    shellString = (
        "python " + run_script + " -y " + str(yco) + " -p " + str(p_bar) + 
        " -t " + str(t_c) + " -s " + str(s_box) + " -m " + str(n_moves) + 
        " -e " + str( n_equil) + " -o " + str(n_prod) + " -f " + str(filepath)
    )
    
    subprocess.run(shellString, shell = True)
    
    return()
    #return( output )

# if __name__ == '__main__':
#     experiments = pd.read_csv(exp_filepath)
    
#     process_list = []
#     for i in range(experiments.shape[0]):
#         p =  Process(target= run_single, args = (experiments,i))
#         p.start()
#         process_list.append(p)
    
#     for process in process_list:
#         process.join()
        
if __name__ == '__main__':
    experiments = pd.read_csv(exp_filepath)
    
    pool = Pool(processes = cpu_count()-1)

    for i in range(experiments.shape[0]):
        pool.apply_async(run_single, args=(experiments,i))
        
    pool.close()
    pool.join()


print( "Done!")
 
# experiments = pd.read_csv(exp_filepath)
       
# for i in range(experiments.shape[0]):
    
#     if experiments.complete[i] == 0:
#         experiments = run_single( experiments, i)
#         experiments.to_csv(exp_filepath, index = False)


