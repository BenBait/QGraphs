cs *multiply_gen (cs *A, cs *B, int n_complete);

// multiply 2 matrices in a dumb way (because it requires checking what index
// we are at). This is a general matrix multiplication
cs *multiply_gen (cs *A, cs *C, int n_complete)
{
    cs *L = malloc(sizeof(cs));
    assert(L != NULL);
    // A and B must be square, just trust without checking for now
    int n = A->n;
    L->m = n;
    L->n = n;
    int num_entries = A->nz;
    L->nz = num_entries;
    L->nzmax = NZMAX;
    L->x = calloc(n*n, sizeof(T));
    L->p = calloc(n*n, sizeof(int));
    L->i = calloc(n*n, sizeof(int));
    assert(L->x != NULL);
    assert(L->p != NULL);
    assert(L->i != NULL);
    printf("A:\n");
    print_dense(A, num_entries, n);
    //printf("B:\n");
    //print_dense(B, num_entries, n);

    // get the transpose of the right matrix, B, then multiply every row down
    // to get the entries in the current row of the left matrix, A
    //cs *C = transpose(B);
    //printf("C:\n");
    //print_dense(C, num_entries, n);

    // id in the output matrix
    int id = 0;
    // each row in the top n_complete rows has n_complete entries, except for
    // the one with the entry we are looking for (which is one less)
    int j = 0;
    int k = 0;
    int r, c;
    int C_r;
    T val = 0;
    int last_col   = -1;
    int top_of_row = 0;
    //for (j = 0; j < num_entries; j++) { // for entries in A
    int rows_of_A = 0;
    while (rows_of_A < n) {
        for (j = 0; j < num_entries; j++) { // for entries in B
            // output matrix entry row is at the row of A
            //printf("Getting R\n");
            r = A->i[k];
            // we've done the transpose of B, so the entry in the output matrix
            // is actually the row of the current version
            //printf("Getting C\n");
            c = C->p[j];
            C_r = C->i[j];
            //printf("Cr: %d\n", C_r);
            //printf("Getting Val\n");
            val = A->x[k] * C->x[j];
            //printf("Got indices\n");
            if (last_col != C_r) { // we've gone to a new row
                // recall id is incremented in push
                //printf("Pushing\n");
                cs_push(L, r, c, val, &id);
                //printf("Pushed\n");
                k = top_of_row;
                last_col = C_r;
            } else {
                // add the value of the next product in the dot product
                // of the rows
                //if (id > num_entries) printf("BAD\n");
                //printf("Incrementing\n");
                L->x[id-1] += val;
                //printf("Incremented\n");
                k++;
            }
        }
        top_of_row = k;
        rows_of_A++;
    }
    /* ATTEMPT 1:
     * indices are not right
    for (i = 0; i < n_complete; i++) {
        // sum up all the products for this entry
        for (j = 0; j < num_row_entries_top; j++) {
            val += A->x[j] * C->x[j];
        }
        assert(A->i[j] == B->i[j]);
        assert(A->p[j] == B->p[j]);
        // output matrix entry row is at the row of A
        r = A->i[j];
        // we've done the transpose of B, so the entry in the output matrix
        // is actually the row of the current version
        c = B->p[j];
        cs_push(L, i, c, val, &i);
    }
    */
}
