/*
 * header file with functions to generate matrices 
 * A is always the adjacency matrix
 */

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

void print_SparseM(SparseM *A, int num_entries, int n) {

    T *full = calloc(n*n, sizeof(T));

    int r, c;
    T val;
    int full_id;

    for (int i = 0; i < num_entries; i++) {
        r = (A->data[i]).r;
        c = (A->data[i]).c;
        val = (A->data[i]).val;
        full_id = c + r*n;
        full[full_id] = val;
    }

    for (int r = 0; r < n; r++) {
        for (int c = 0; c < n; c++) {
            full_id = c + r*n;
            printf("%.2f\t", full[full_id]);
        }
        printf("\n");
    }
}

SparseM *gen_from_text(char *filename, int *num_entries, int *n) {

    FILE *infile = fopen(filename, "r");
    *n           = fgetc(infile);
    *num_entries = fgetc(infile);

    SparseM *A = malloc(sizeof(SparseM));
    assert(A != NULL);
    // allocate the data in the matrix
    A->data = calloc(*num_entries, sizeof(SparseM_entry));
    assert(A->data != NULL);

    int i = 0;
    for (i = 0; i < *num_entries; i++) {
        (A->data[i]).r = fgetc(infile);
        (A->data[i]).c = fgetc(infile);
        (A->data[i]).val = fgetc(infile);
    }

    fclose(infile);
    return A;
}

// matrix to generate the identity matrix in sparse format
void get_identity_matrix(SparseM *A, int n) {

    // make identity matrix
    for (int i = 0; i < n; i++) {
        // NEED BOUND CHECKING
        (A->data[i]).r = i;
        (A->data[i]).c = i;
        (A->data[i]).val = 1;
    }
}

// matrix to generate the identity matrix in sparse format
SparseM *get_fully_connected(int num_entries, int n) {

    SparseM *A = malloc(sizeof(SparseM));
    assert(A != NULL);
    // allocate the data in the matrix
    A->data = calloc(num_entries, sizeof(SparseM_entry));
    assert(A->data != NULL);

    // better way to do this with decreasing index
    // this will generate the adjacency matrix of a fully connected graph
    int idx = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // NEED BOUND CHECKING
            if (i != j) {
                (A->data[idx]).r = i;
                (A->data[idx]).c = j;
                (A->data[idx]).val = 1;
                idx++;
            }
        }
    }
    assert(idx == num_entries);

    return A;
}

// idx keeps track of the number of entries in the sparse matrix
void add_path(SparseM *A, int n, int n_complete, int start_N, int end_N, int *next_N, T p, int *idx) {

    printf("ADDING PATH: n is %d\n", n);
    printf("NEXT: %d\n", *next_N);
    // connect the start node with the next node, which is the first of the 2
    // nodes in the segmented edge
    int r = start_N;
    int c = *next_N;
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
    // we might increment to the 0 node again
    c = (r+1) % n;
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
    (A->data[idy]).r = r;
    (A->data[idy]).c = c;
    (A->data[idy]).val = 1 - p;
    idy++;

    // now do the reverse direction for the connection to the perimeter
    r = end_N;
    c = ((*next_N)+1) % n;
    (A->data[idy]).r = r;
    (A->data[idy]).c = c;
    (A->data[idy]).val = 1 / ((float)n_complete-1);
    idy++;

    r = *next_N;
    c = start_N;
    (A->data[idy]).r   = r;
    (A->data[idy]).c   = c;
    // chance of returning is 1 - p
    (A->data[idy]).val = 1 - p;
    idy++;

    // get the updated current index
    *idx = idy;
    // increment next_N by 3 for next time
    *next_N += 3;
    print_SparseM(A, 54, n);
}

SparseM *get_segmented_complete(int num_entries, int n, T p, int segs) {

    printf("NUM ENTRIES: %d\n", num_entries);
    SparseM *A = malloc(sizeof(SparseM));
    assert(A != NULL);
    // allocate the data in the matrix
    A->data = calloc(num_entries, sizeof(SparseM_entry));
    assert(A->data != NULL);

    int j,k;
    int n_complete = n / segs;
    // number of paths between an orig node and the other orig nodes
    // that are not the adjacent original nodes, these are "inside" paths
    int n_in = n_complete - 3;
    int end_N;

    // next_N is the index of the first node in the segmented edge after the
    // orig node
    int next_N = 1;

    // idx keeps track of the number of entries in the sparse matrix
    int idx = 0;

    // add the segmented edges between all of the perimeter nodes
    for (j = 0; j < n; j += segs) {
        add_path(A, n, n_complete, j, (j+segs) % n, &next_N, p, &idx);
    }

    //printf("PRINTING IN GENERATION FUNCTION\n");
    //print_SparseM(A, num_entries, n);

    // The outside nodes are multiples of 3
    // Add the paths between the first 2 non sequential perimeter nodes
    for (j = 0; j < 6; j+=segs) {
        // j is the start node
        for (k = 2; k < n_in+2; k++) {
            // k calculates the end node, and it starts at 2, so it ends at the
            // number of inside paths + 2
            end_N = k*segs + j;
            add_path(A, n, n_complete, j, end_N, &next_N, p, &idx);
        }
    }

    // increment this to avoid duplicate paths
    int m = 1;
    for (j = 6; j < n; j += segs) {
        // j is the start node
        for (k = 2; k < (n_in-m)+2; k++) {
            // k calculates the end node, and it starts at 2, so it ends at the
            // number of inside paths + 2
            end_N = k*segs + j;
            add_path(A, n, n_complete ,j, end_N, &next_N, p, &idx);
        }
        m++;
    }

    return A;
}

//SOME USAGE NOTES
    /*
    if (argc == 2) {
        SparseM *A = gen_from_text(argv[1], &num_entries, &n);
    } else {
        SparseM *A = get_fully_connected(num_entries, n);
    }*/
