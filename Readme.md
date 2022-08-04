# Get Graph Laplacian

FILES: get\_laplacian.c -> converts adj and deg matrices to laplacian; create\_matrix.c -> creates input adj and deg matrices

## Command Line Input

### Input Matrix Format:
 - {row}{column}{value}
 - Each row and column field in {} is B bits where B = ceil (log<sub>2<\sub>(N))
 - Number of bits for value is defined in header (Everything else is zero)
 - Most nodes nodes are not connected, so we use a sparse matrix format

### Header:
 - {N}, {Number of Value Bits}

### Input Background Info:
 - Adjacency matrix representation of a graph
 - Matrix with N nodes ==> NxN matrix whose entries are u<sub>rc</sub> and
u<sub>rc</sub> = 1 if there is a path between node r and node c

## TO DO:

    1. Separate file to generate matrices and output sparse matrix

    2. I/O: getLaplacian executable gets sparse matrix

    3. I/O: getLaplacian executable prints sparse matrix

    4. Function: input matrices (adjacency and degree) --> graph Laplacian
