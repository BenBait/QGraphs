# Get Graph Laplacian

build with `gcc get_laplacian.c`  

## Command Line Input
Adjacency matrix representation of a graph

Matrix with N nodes ==> NxN matrix whose entries are u<sub>rc</sub> and
u<sub>rc</sub> = 1 if there is a path between node r and node c, the entry is 1

Most nodes nodes are not connected, hence we will use a sparse matrix format

Input Matrix Format: {row}{column}{value}
Each row and column field in {} is B bits where B = ceil (log<sub>2<\sub>(N))
Number of bits for value is defined in header
(Everything else is zero)

Header: 2 integers. {N}, {Number of Value Bits}


TO DO:

    1. Separate file to generate matrices and output sparse matrix

    2. I/O: get Laplacian executable gets sparse matrix

    3.  
