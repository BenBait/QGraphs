# data and visualization
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

# scripting capabilities
import os

command = "./build.sh 2"
os.system(command)

# number of nodes in the complete graph to try
num_nodes = [3, 4, 5]
# possible transition probabilities, p, to use
ps        = [0.5, 0.2, 0.333]
# possible target nodes, for now just the 1st
targets   = [1]
times  = [1, 10, 100, 300]
gammas = [1, 1.5, 2, 2.09]
#call the data gathering C executable while iterating through possible params
for n in num_nodes:
    for p in ps:
        for w in targets:
            for t in times:
                for g in gammas:
                    params = f"{n} {p} {w} {t} {g}"
                    command = f"./a.out {params}"
                    os.system(command)
