import subprocess
import pandas as pd
#import numpy as np
#import sys
import os
from time import time
from multiprocessing import Pool, cpu_count


exp_filepath = "carbon_wall/distribution_test/distribution_test_exp.csv"
data_filepath = "carbon_wall/distribution_test/gas"
run_script = "carbon_wall/distribution_test/distribution_runs.py"

# print( "Available CPU:", cpu_count() )

def run_single( exp_list, i):
    yco = exp_list.yco[i]
    p_bar = exp_list.p_bar[i] 
    t_c = exp_list.t_c[i]
    s_box = exp_list.s_box[i]
    n_moves = int( exp_list.n_moves[i] )
    n_equil = int( exp_list.n_equil[i] )
    n_prod = int( exp_list.n_prod[i] )
    filepath = data_filepath + "/" + str(exp_list.exp[i])
    sf = exp_list.sf[i]
    W = exp_list.W[i]
    #e_co = int( exp_list.e_co[i] )
    #s_co = int( exp_list.s_co[i] )
    #e_me = int( exp_list.e_me[i] )
    #s_me = int( exp_list.s_me[i] )
    
    shellString = (
        "python " + run_script + " -y " + str(yco) + " -p " + str(p_bar) + 
        " -t " + str(t_c) + " -s " + str(s_box) + " -m " + str(n_moves) + 
        " -e " + str( n_equil) + " -o " + str(n_prod) + " -f " + str(filepath)  +
        # " --e_co " + str(e_co) + " --s_co " + str(s_co) +
        # " --e_me " + str(e_me) + " --s_me " + str(s_me) 
        " --sf " + str(sf) + " --w " + str(W)
    )
    
    t0 = time()
    processOutput = subprocess.run(shellString, shell = True, text=True, stdout=subprocess.PIPE)
    tf = time()
    
    runtime = round( (tf - t0)/60, 2 )
    
    return( "Run " + str( i ) + " complete in " + str(runtime) + " mins!" )
    #return( output )
    
    
def print_output(result):
    print( result )
        
if __name__ == '__main__':
    experiments = pd.read_csv(exp_filepath)
    
    pool = Pool(processes = cpu_count()-1)

    for i in range(experiments.shape[0]):
        pool.apply_async(run_single, args=(experiments,i),callback=print_output)
        
    pool.close()
    pool.join()

#if __name__ == '__main__':
#    experiments = pd.read_csv(exp_filepath)
#    runOutput = run_single( experiments,  4 )
#    print( runOutput.stdout )



print( "Done!")
 
# experiments = pd.read_csv(exp_filepath)
       
# for i in range(experiments.shape[0]):
    
#     if experiments.complete[i] == 0:
#         experiments = run_single( experiments, i)
#         experiments.to_csv(exp_filepath, index = False)


