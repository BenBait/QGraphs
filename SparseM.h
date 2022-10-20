#include<math.h>
/*
 * header file with functions to generate matrices 
 * A is always the adjacency matrix
 */

// floating point epsilon. use 10^-15
#define EPSILON 0.0000000000000001
// data type for matrix entries
typedef float T;

// entry into a sparse matrix (only used for nonzero values)
typedef struct SparseM_entry {
    // row and column of nonzero entry
    int r,c;
    // nonzero entry
    T val;
} SparseM_entry;

// array of SparseM_entries with nonzero data
typedef struct SparseM {
    int num_entries;
    SparseM_entry *data;
} SparseM;

void print_SparseM (SparseM *A, int num_entries, unsigned long n)
{

    T *full = calloc(n*n, sizeof(T));

    int r, c;
    T val;
    int full_id;

    for (int i = 0; i < num_entries; i++) {
        r = (A->data[i]).r;
        c = (A->data[i]).c;
        val = (A->data[i]).val;
        full_id = c + r*n;
        assert(full_id < n*n);
        full[full_id] = val;
    }

    for (int r = 0; r < n; r++) {
        for (int c = 0; c < n; c++) {
            full_id = c + r*n;
            if (fabs(full[full_id]) < EPSILON)
                printf("%d\t", 0);
            else
                printf("%.2f\t", full[full_id]);
        }
        printf("\n");
    }

    free(full);
}

/* sum all the rows of a sparse matrix and confirm the each sum to one
 * this test does not scale well for large n
 */
void dumb_row_sum_SparseM (SparseM *A, int num_entries, unsigned long n, int P)
{
    T *full = calloc(n*n, sizeof(T));

    int r, c;
    T val;
    unsigned long full_id;

    // n^2 might be too big for a regular int (if n_complete>=220 this happens)
    unsigned long n_sq = (unsigned long) n*n;

    for (int i = 0; i < num_entries; i++) {
        r = (A->data[i]).r;
        c = (A->data[i]).c;
        val = (A->data[i]).val;
        full_id = c + r*n;

        assert(full_id < n_sq);
        full[full_id] = val;
    }

    float row_sum;
    for (int r = 0; r < n; r++) {
        row_sum = 0;
        for (int c = 0; c < n; c++) {
            full_id = c + r*n;
            row_sum += full[full_id];
        }
        // assert sum is within floating point error of 1
        if (P) printf("row %d sums to %.1f...\n", r, row_sum);
        assert(abs(row_sum - 1.0) < EPSILON);
    }

    free(full);
}

// idx keeps track of the number of entries in the sparse matrix
void add_path (SparseM *A, unsigned long n, int n_complete, int num_entries, int start_N, int end_N, int next_N, T p, int *idx)
{
    // connect the start node with the next node, which is the first of the 2
    // nodes in the segmented edge
    int r = start_N;
    int c = next_N;
    assert(r < n);
    assert(c < n);
    // idy stores the current id
    int idy = *idx;
    (A->data[idy]).r   = r;
    (A->data[idy]).c   = c;
    // adjacent node to perimeter node has equal chance of getting traversed
    // as the other nodes
    (A->data[idy]).val = 1 / ((float)n_complete-1);
    idy++;

    // connect first node of the edge with the second
    r = c;
    c = (r+1) % n;
    assert(r < n);
    assert(c < n);
    (A->data[idy]).r = r;
    (A->data[idy]).c = c;
    (A->data[idy]).val = p;
    idy++;
    // these nodes have the same probability in both directions
    (A->data[idy]).r = c;
    (A->data[idy]).c = r;
    (A->data[idy]).val = p;
    idy++;

    // connect second node on the edge with the terminating perimeter node
    r = c;
    c = end_N;
    assert(r < n);
    assert(c < n);
    (A->data[idy]).r = r;
    (A->data[idy]).c = c;
    (A->data[idy]).val = 1 - p;
    idy++;
 
    // now do the reverse direction for the connection to the perimeter
    r = end_N;
    c = (next_N+1) % n;

    assert(r < n);
    assert(c < n);
    (A->data[idy]).r = r;
    (A->data[idy]).c = c;
    (A->data[idy]).val = 1 / ((float)n_complete-1);
    idy++;

    r = next_N;
    c = start_N;
    assert(r < n);
    assert(c < n);
    (A->data[idy]).r   = r;
    (A->data[idy]).c   = c;
    // chance of returning is 1 - p
    (A->data[idy]).val = 1 - p;
    idy++;

    // get the updated current index
    *idx = idy;
}

/* generates the adjacency matrix for a complete matrix with n_complete entries
 * and is then segmented */
SparseM *get_segmented_complete (int num_entries, int n_complete, T p, int segs)
{
    SparseM *A = malloc(sizeof(SparseM));
    assert(A != NULL);
    // allocate the data in the matrix
    A->data = calloc(num_entries, sizeof(SparseM_entry));
    assert(A->data != NULL);

    int j,k;
    // number of paths between an orig node and the other orig nodes
    // that are not the adjacent original nodes, these are "inside" paths
    int n_in = n_complete - 3;
    int end_N;

    // next_N is the index of the first node in the segmented edge after the
    // orig node
    int next_N = 1;

    // idx keeps track of the number of entries in the sparse matrix
    int idx = 0;

    int num_paths = 0;
    // add the segmented edges between all of the original nodes
    for (j = 0; j <= (n_complete-1)*segs; j += segs) {
        add_path(A, n, n_complete, num_entries, j, (j+segs) % (n_complete*segs), next_N, p, &idx);
        next_N += 3;
        num_paths++;
    }
    // the next node is the first in the 0 <-> 6 path (for all n), so we went past it
    next_N -= 1;
    /*
    // Optional Test part 1
    assert(num_paths == n_complete);
    */

    // The outside nodes are multiples of 3
    // Add the paths between the first 2 non sequential perimeter nodes
    for (j = 0; j < 2*segs; j+=segs) {
        // j is the start node
        for (k = 2; k < n_in+2; k++) {
            // k calculates the end node, and it starts at 2, so it ends at the
            // number of inside paths + 2
            end_N = k*segs + j;
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
    for (j = segs*2; j <= (n_complete-1)*segs; j += segs) {
        // j is the start node
        for (k = 2; k < (n_in-m)+2; k++) {
            // k calculates the end node, and it starts at 2, so it ends at the
            // number of inside paths + 2
            end_N = k*segs + j;
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
int *get_degree(SparseM *A, int num_entries, int n) {

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

    for (i = 0; i < num_entries; i++) {
        // get nonzero entries
        r = (A->data[i]).r;
        // don't actually need the column here

        // accumulate the degree for each entry
        if (A->data[i].val > EPSILON) {
            d[r] += 1;
        }
    }

    return d;
}

/* Laplacian is the degree matrix minus the adjacency matrix */
SparseM *get_laplacian(SparseM *A,
                       int num_entries,
                       int n) {

    // this is the Laplacian we are going to use
    SparseM *L = malloc(sizeof(SparseM));
    assert(L != NULL);

    // entries in the laplacian are the entries in A, plus the diagonals
    int num_L_entries = num_entries + n;

    L->data = calloc(num_L_entries, sizeof(SparseM_entry));
    assert(L->data != NULL);

    // loop counter
    int i = 0;
    
    // get the diagonals of the degree matrix
    int *d = get_degree(A, num_entries, n);
    
    // row and column indices
    int r = 0;
    int c = 0;

    // index in the adjacency matrix
    int idx = 0;

    // L = D - A, D is diagonal
    // diagonals of A are always 0, so just fill those with the values of d
    for (i = 0; i < n; i++) {
        // note that this only iterates through n of num_L_entries
        (L->data[i]).r = i;
        (L->data[i]).c = i;
        (L->data[i]).val = d[i];
    }
 
    // NEED TO ADD BOUNDS CHECKING: NOT MORE THAN N VALUES FOR ANY GIVEN R
    for (i = n; i < num_L_entries; i++) {
        // get nonzero entries
        // indices in A are from 0 to num_entries
        (L->data[i]).r = (A->data[i-n]).r;
        (L->data[i]).c = (A->data[i-n]).c;
        (L->data[i]).val = -((A->data[i-n]).val);
    }

    free(d);
    return L;
}
