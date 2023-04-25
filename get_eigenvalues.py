import math
import numpy as np
from numpy.linalg import eig
import os

#A = [[1, -0.3], [-0.3, 1]]
#w,v = eig(A)
#print(f"eigenvalues: {w}")
#print(f"eigenvectors: {v}")
#exit()

cmd = "./build.sh 2"
os.system(cmd)

#num_nodes = [4, 5, 6]
num_nodes = [3, 4, 5, 6, 7]
#ps        = [0.1, 0.3, 0.5]
ps        = [0.1]

def kth_eigenvector(n, k):
    if (k < 1):
        print("k must be greater than 0")
        exit()
    v = np.zeros(n)
    v[k]   =  1 / math.sqrt(2)
    v[k-1] = -1 / math.sqrt(2)
    return v

for n in num_nodes:
    for p in ps:
        cmd = f"./a.out {n} {p} 1 1 1"
        os.system(cmd)

        L = np.genfromtxt("laplacian.matrix", delimiter="\t", dtype=np.double)
        #print(L)
        N = L.shape[0]
        print(N)
        #A = L[np.ix_([0,n],[0,n])]
        A = L[np.ix_(np.arange(0,n),np.arange(0,n))]
        '''
        B = L[np.ix_(np.arange(0,n),np.arange(n,N-n))]
        C = L[np.ix_(np.arange(n,N-n),np.arange(0,n))]
        D = L[np.ix_(np.arange(n,N-n),np.arange(n,N-n))]
        '''
        B = L[np.ix_(np.arange(0,n),np.arange(n,N))]
        C = L[np.ix_(np.arange(n,N),np.arange(0,n))]
        D = L[np.ix_(np.arange(n,N),np.arange(n,N))]
        print(B)
        continue

        w,v = eig(L)

        dbl_prec = np.finfo(np.double).precision
        epsilon = 10**(-dbl_prec)

        # round down to zero for small values
        w[abs(w - 0) < epsilon] = 0
        w.sort()
        np.set_printoptions(suppress=True)
        print(f"eigenvalues: {w}\n")
        #print(f"eigenvectors: {v}")
        
        Z = np.identity(N-n) * w[0]
        diff = D - Z
        neg_inv_diff = -1*np.linalg.inv(diff)

        v_0 = kth_eigenvector(n, 1)
        Cv_0 = C.dot(v_0)

        v_prime = neg_inv_diff.dot(Cv_0)
        print("v_prime:")
        print(v_prime)

        v = np.append(v_0, v_prime)

        print("L*v:")
        print(L.dot(v))
