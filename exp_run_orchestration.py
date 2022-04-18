import subprocess
import pandas as pd
#import numpy as np
#import sys
#import os
from multiprocessing import Pool, cpu_count


exp_filepath = "longer_tests/c02_pure_component_exp_list.csv"
data_filepath = "longer_tests/"
run_script = ""

# print( "Available CPU:", cpu_count() )

def run_single( exp_list, i):
    yco = exp_list.yco[i]
    p_mpa = exp_list.p_mpa[i] 
    t_c = exp_list.t_c[i]
    s_box = exp_list.s_box[i]
    n_moves = exp_list.n_moves[i]
    n_equil = exp_list.n_equil[i]
    n_prod = exp_list.n_prod[i]
    filepath = data_filepath + str(exp_list.exp[i]) + ".csv"
    
    shellString = (
        "python " + run_script + " -y " + str(yco) + " -p " + str(p_mpa) + 
        " -t " + str(t_c) + " -s " + str(s_box) + " -m " + str(n_moves) + 
        " -e " + str(n_equil) + " -o " + str(n_prod) + " -f " + str(filepath)
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
    
    pool = Pool(processes = cpu_count()-2)

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


