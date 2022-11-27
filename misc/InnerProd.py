import numpy as np
from scipy import linalg


f = open("tri.txt", "r")
L = np.genfromtxt(f, delimiter=",")
print(L)

f1 = open("m3.txt", "r")
M = np.genfromtxt(f1, delimiter=",", dtype=np.double)
print(M)

f2 = open("m3_inv.txt", "r")
M_inv = np.genfromtxt(f2, delimiter=",", dtype=np.double)
print(M_inv)

L_dag_tmp = L.dot(M_inv)
print(L_dag_tmp)
L_dag = M.dot(L_dag_tmp)
print(L_dag)
