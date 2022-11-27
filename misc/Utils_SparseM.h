/* Misc functions that can be utilized in SparseM.h */

SparseM *gen_from_text (char *filename, int *num_entries, int *n) 
{

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
void get_identity_matrix (SparseM *A, unsigned long n)
{

    // make identity matrix
    for (int i = 0; i < n; i++) {
        // NEED BOUND CHECKING
        (A->data[i]).r = i;
        (A->data[i]).c = i;
        (A->data[i]).val = 1;
    }
}

// matrix to generate the identity matrix in sparse format
SparseM *get_fully_connected (int num_entries, unsigned long n)
{

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

