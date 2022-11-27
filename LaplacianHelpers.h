#ifndef L_HELP
#define L_HELP
/*
 * header file with functions to generate matrices 
 * A is always the adjacency matrix
 */
#include<math.h>
#include"SparseMethods.h"


// idx keeps track of the number of entries in the sparse matrix
//void add_path (cs *A, unsigned long n, int n_complete, int num_entries,
               //int start_N, int end_N, int next_N, T p)
void add_path (cs *A, unsigned long n, int n_complete, int num_entries,
               int start_N, int end_N, int next_N, T p, int *idx)
{
    // connect the start node with the next node, which is the first of the 2
    // nodes in the segmented edge
    int r = start_N;
    int c = next_N;
    assert(r < n && c < n);

    int idy = *idx;
    cs_push(A, r, c, 1 / ((T)n_complete-1), &idy);
    // connect first node of the edge with the second

    r = c;
    c = (r+1) % n;
    assert(r < n && c < n);
    cs_push(A, r, c, p, &idy);
    // these nodes have the same probability in both directions
    cs_push(A, c, r, p, &idy);

    // connect second node on the edge with the terminating perimeter node
    r = c;
    c = end_N;
    assert(r < n && c < n);
    cs_push(A, r, c, 1-p, &idy);
 
    // now do the reverse direction for the connection to the perimeter
    r = end_N;
    c = (next_N+1) % n;
    assert(r < n && c < n);
    cs_push(A, r, c, 1 / ((T)n_complete-1), &idy);

    r = next_N;
    c = start_N;
    assert(r < n && c < n);
    cs_push(A, r, c, 1-p, &idy);

    // get the updated current index
    *idx = idy;
}

/* generates the adjacency matrix for a complete matrix with n_complete entries
 * and is then segmented */
cs *get_segmented_complete (int num_entries, int n, int n_complete, T p, int segs)
{
    cs *A = malloc(sizeof(cs));
    assert(A != NULL);
    // square matrix; m = n = number of columns
    A->m = n;
    A->n = n;
    A->nz = num_entries;
    A->nzmax = NZMAX;
    // allocate the data in the matrix
    A->x = calloc(num_entries, sizeof(T));
    // column indices
    A->p = calloc(num_entries, sizeof(int));
    // row indices
    A->i = calloc(num_entries, sizeof(int));
    assert(A->x != NULL);
    assert(A->p != NULL);
    assert(A->i != NULL);

    int j,k;
    // number of paths between an orig node and the other orig nodes
    // that are not the adjacent original nodes, these are "inside" paths
    int n_in = n_complete - 3;
    int end_N;

    // next_N is the index of the first node after the n_complete nodes
    int next_N = n_complete;

    // idx keeps track of the number of entries in the sparse matrix
    int idx = 0;

    int num_paths = 0;
    // add the segmented edges between all of the original nodes
    for (j = 0; j < n_complete; j++) {
        // for 3 segs per edge, a path will first be added between node 0 and node 1
        add_path(A, n, n_complete, num_entries, j, (j+1) % n_complete, next_N, p, &idx);
        next_N += 2;
        num_paths++;
    }
    /*
    // Optional Test part 1
    assert(num_paths == n_complete);
    */

    // The outside nodes are multiples of 3
    // Add the paths between the first 2 non sequential perimeter nodes
    for (j = 0; j < 2; j++) {
        // j is the start node
        for (k = 2; k < n_in+2; k++) {
            // k calculates the end node, and it starts at 2, so it ends at the
            // number of inside paths + 2
            end_N = k + j;
            add_path(A, n, n_complete, num_entries, j, end_N, next_N, p, &idx);
            num_paths++;
            next_N += 2;
        }
    }
    /*
    // Optional Test part 2
    // number of paths through the inside, for now just the first 2 nodes
    int num_inside = 2*(n_complete - segs);
    assert(num_paths == (n_complete + num_inside));
    */

    // increment this to avoid duplicate paths
    int m = 1;
    for (j = 2; j < n_complete; j++) {
        // j is the start node
        for (k = 2; k < (n_in-m)+2; k++) {
            // k calculates the end node, and it starts at 2, so it ends at the
            // number of inside paths + 2
            end_N = k + j;
            // the final "next N" will be n, which is too big
            add_path(A, n, n_complete, num_entries, j, end_N, next_N, p, &idx);
            next_N += 2;
            num_paths++;
        }
        m++;
    }
    /*
    // Optional Test part 3
    // now count the remaining inside edges
    int added_edges = n_complete - segs - 1;
    while (added_edges > 0) {
        num_inside += added_edges;
        added_edges--;
    }
    assert(num_paths == (n_complete + num_inside));
    */

    return A;
}

/* Helper function for getting the laplacian */
int *get_degree(cs *A, int num_entries, int n) {

    // loop counter
    int i = 0;
    // row and column counters
    int r = 0;
    int c = 0;
    // d contains the diagonal elements of the degree matrix
    // i.e. d[i] = D_ii
    int *d = calloc(n, sizeof(T));

    // accumulate the value of d for each row in this var
    int curr_d = 0;
    int idx = 0;

    for (i = 0; i < A->nz; i++) {
        // get nonzero entries
        r = A->i[i];
        // don't actually need the column here

        // accumulate the degree for each entry
        if (creal(A->x[i]) > 0) {
            d[r] += 1;
        }
    }

    return d;
}

/* Laplacian is the degree matrix minus the adjacency matrix */
cs *get_laplacian(cs *A, int num_entries, int n)
{
    // this is the Laplacian we are going to use
    cs *L = malloc(sizeof(cs));
    assert(L != NULL);

    // entries in the laplacian are the entries in A, plus the diagonals
    int num_L_entries = num_entries + n;
    // square matrix; m = n = number of columns
    L->m = n;
    L->n = n;
    L->nz = num_L_entries;
    L->nzmax = NZMAX;
    // allocate the data in the matrix
    L->x = calloc(num_L_entries, sizeof(T));
    // column indices
    L->p = calloc(num_L_entries, sizeof(int));
    // row indices
    L->i = calloc(num_L_entries, sizeof(int));
    assert(L->x != NULL);
    assert(L->p != NULL);
    assert(L->i != NULL);

    // get the diagonals of the degree matrix
    int *d = get_degree(A, num_entries, n);

    // row and column indices
    int r = 0;
    int c = 0;
    T val = 0;

    // L = D - A, D is diagonal
    // diagonals of A are always 0, so just fill those with the values of d
    int i = 0;
    while (i < n) {
        // note that this only iterates through n of num_L_entries
        
        // i is iterated in cs_push
        // use actual d for sanity checking the sum of each row is 0
        //cs_push(L, i, i, (T)d[i], &i);
        cs_push(L, i, i, 1.0, &i);
    }
 
    // NEED TO ADD BOUNDS CHECKING: NOT MORE THAN N VALUES FOR ANY GIVEN R
    while (i < num_L_entries) {
        // get nonzero entries
        // indices in A are from 0 to num_entries
        r = A->i[i-n];
        c = A->p[i-n];
        val = -(A->x[i-n]);
        cs_push(L, r, c, val, &i);
    }

    free(d);
    return L;
}

/* Laplacian is the degree matrix minus the adjacency matrix */
cs *get_hamiltonian (cs *A, int n, int w, T gamma)
{
    // this is the Laplacian we are going to use
    cs *L = malloc(sizeof(cs));
    assert(L != NULL);

    // entries in the laplacian are the entries in A, plus the diagonals
    int num_L_entries = A->nz + n;
    // square matrix; m = n = number of columns
    L->m = n;
    L->n = n;
    L->nz = num_L_entries;
    L->nzmax = NZMAX;
    // allocate the data in the matrix
    L->x = calloc(num_L_entries, sizeof(T));
    // column indices
    L->p = calloc(num_L_entries, sizeof(int));
    // row indices
    L->i = calloc(num_L_entries, sizeof(int));
    assert(L->x != NULL);
    assert(L->p != NULL);
    assert(L->i != NULL);

    // get the diagonals of the degree matrix
    int *d = get_degree(A, A->nz, n);

    // row and column indices
    int r = 0;
    int c = 0;
    T val = 0;

    // L = D - A, D is diagonal
    // diagonals of A are always 0, so just fill those with the values of d
    // to sanity checking the sum of each row is 0 : cs_push(L,i,i,(T)d[i],&i);
    // BUT we want the probabilistic laplacian, so the diagonals are just 1

    // note that this only iterates through n of num_L_entries
    int i = 0;
    while (i < n) {
        // i is iterated in cs_push

        // to get the Hamiltonian, multiply the Laplacian by gamma,
        // then sub the oracle function
        val = gamma;
   
        // if we are at the search node, subtract the value in the oracle
        // function from the scaled laplacian
        if (i == w) val -= 1.0;

        cs_push(L, i, i, val, &i);
    }
 
    // NEED TO ADD BOUNDS CHECKING: NOT MORE THAN N VALUES FOR ANY GIVEN R
    while (i < num_L_entries) {
        // get nonzero entries
        // indices in A are from 0 to num_entries
        r = A->i[i-n];
        c = A->p[i-n];
        val = gamma * -(A->x[i-n]);

        if (i == w) val -= 1.0;

        cs_push(L, r, c, val, &i);
    }

    free(d);
    return L;
}

#endif
