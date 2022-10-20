import argparse
import pandas as pd
import numpy as np
exec(open("gcmc.py").read())

yco = 0
P_bar = 10 #bar
T_C = 45 #bar
s_box = 34 #A
n_moves = 1000
n_equil = 100000
n_prod = 1