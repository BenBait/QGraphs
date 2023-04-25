# This python program creates the Laplacian for the segmented-complete graph
# See section 3
import numpy as np

# helper function for calculating the B submatrix for the laplacian
def write_circulent_submatrices(N0: int, N: int, q: float):
    B = np.zeros((N0, N-N0))
    row = np.zeros(N-N0)
    row[0] = row[-1] = q
    B[0] = row
    for i in range(1, N0):
        row = np.roll(row, 1)
        B[i] = row
    return B

import numpy as np

import numpy as np

def generate_row(N, i):
    print(N)
    row = np.zeros(N**2 - N)
    k = i * (N-1) * 2
    for j in range(N-2):
        row[(k+j)] = -1 / (N-1)
        row[(k+j+(N-1)*2)] = -1 / (N-1)
    return row

def create_B(N):
    # Initialize an N x (2N-1) matrix of zeros
    matrix = np.zeros((N, N**2 - N))
    
    # Populate the matrix
    for i in range(N):
        matrix[i] = generate_row(N, i)

    return matrix

def get_segmented_complete_laplacian(N0: int, p: float):
    # segments on each edge of the complete graph
    segs = 3;
    new_nodes_per_edge = 2; # <== there being 3 segs per node

    # number of edges in the complete graph
    E0 = int(N0*(N0 - 1) / 2)

    # number of nodes in the segmented complete graph
    N = N0 + 2*E0

    # Create A matrix (N0 x N0 identity)
    A = np.identity(N0)

    # Set q value for B
    q = -1/(N0-1)

    # Create B matrix
    q = -1 / (N0-1)
    #B = write_circulent_submatrices(N0, N, q)
    B = create_B(N0)
    print(B)
    # Create C matrix
    q = p - 1
    C = (write_circulent_submatrices(N0, N, q)).T
    print(C)
    exit()
    
    D = np.zeros((N-N0, N-N0))
    for i in range((N-N0)//2):
        D[2*i:2*i+2, 2*i:2*i+2] = np.array([[1, -p], [-p, 1]])
    
    # Create L matrix
    L = np.zeros((N,N))
    L[:N0,:N0] = A
    L[:N0,N0:] = B
    L[N0:,:N0] = C
    L[N0:,N0:] = D
    
    return L

if __name__ == '__main__':
    print(get_segmented_complete_laplacian(4, 0.5))
