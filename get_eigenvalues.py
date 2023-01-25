import numpy as np
from numpy.linalg import eig

#A = [[1, -0.3], [-0.3, 1]]
#w,v = eig(A)
#print(f"eigenvalues: {w}")
#print(f"eigenvectors: {v}")
#exit()

L = np.genfromtxt("laplacian.matrix", delimiter="\t", dtype=np.double)
print(L)

w,v = eig(L)
print(f"eigenvalues: {w}")
print(f"eigenvectors: {v}")
