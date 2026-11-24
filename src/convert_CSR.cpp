// Read .mtx input file and convert to CSR format

#include "SpMM.hpp"

struct edgeData {
    int from;
    int to;
    int val;
};
typedef struct edgeData edgeData;


struct edgeData *queues[10]; 
int queue_Offsets[10];

void make_CSR(const char *input_File);
void make_synthetic(const char *input_File);
void sort_From(edgeData *edges, int nr, int ne);
void sort_To(edgeData *edges, int nc, int ne);

void convert_CSR(const char* name) {
    char input[64];
    strcpy(input, name);
    strcat(input, ".csr");
    
    // if (0) {
    if (access(input, F_OK) != -1) {
        FILE *fp = fopen(input, "rb");
        fread(&mat.nr, sizeof(DATA_TYPE), 1, fp);
        fread(&mat.nc, sizeof(DATA_TYPE), 1, fp);
        fread(&mat.ne, sizeof(DATA_TYPE), 1, fp);
        
        mat.nnz = (DATA_TYPE*)malloc(mat.nr * sizeof(DATA_TYPE));
        mat.ptr = (DATA_TYPE*)malloc((mat.nr + 1) * sizeof(DATA_TYPE));
        mat.idx = (DATA_TYPE*)malloc(mat.ne * sizeof(DATA_TYPE));
        fread(mat.nnz, sizeof(DATA_TYPE), mat.nr, fp);
        fread(mat.ptr, sizeof(DATA_TYPE), mat.nr + 1, fp);
        fread(mat.idx, sizeof(DATA_TYPE), mat.ne, fp);

        fclose(fp);
        
        puts("\n --- Input Matrix Information -------------------------------------");
        printf("| Number of Rows:\t%d\n| Number of Edges:\t%d\n\n", mat.nr - 1, mat.ne);
    }
    else {
        make_CSR(name);
    }
}


void make_CSR(const char *input_File) {
    FILE *fp = fopen(input_File, "rt");

    char line[200];
    bool symmetric = false;
    if (strstr(input_File, ".dat") == NULL) {
        fgets(line, 200, fp);
        if(strstr(line, "symmetric")) symmetric = true;
    }
    //else puts("Synthetic");

    while (fgets(line, 200, fp) != NULL) {
        if (line[0] == '%' || line[0] == '#') continue;   // Comment
        char delims[] = "\t, ";
        char *t;
        t = strtok(line, delims);
        if (t == NULL) {
            fprintf(stderr, "Wrong Format.\n");
            exit(-1);
        }
        mat.nr = atoi(t);
        
        t = strtok(NULL, delims);
        if (t == NULL) {
            fprintf(stderr, "Wrong Format.\n");
            exit(-1);
        }
        mat.nc = atoi(t);
        // assert(mat.nr == atoi(t) + 1);

        t = strtok(NULL, delims);
        if (t == NULL) {
            fprintf(stderr, "Wrong Format.\n");
            exit(-1);
        }
        mat.ne = atoi(t);
        if (symmetric == true) mat.ne *= 2;

        break;
    }

    mat.nnz = (DATA_TYPE*)malloc((mat.nr + 1) * sizeof(DATA_TYPE));
    mat.ptr = (DATA_TYPE*)malloc((mat.nr + 1) * sizeof(DATA_TYPE));
    mat.idx = (DATA_TYPE*)malloc(mat.ne * sizeof(DATA_TYPE));

    edgeData *edges = (edgeData*)malloc(sizeof(edgeData) * mat.ne);
    for (int q = 0; q < 10; ++q) queues[q] = (edgeData*)malloc(sizeof(edgeData) * mat.ne);
    int offset = 0;

    while (fgets(line, 200, fp) != NULL) {
        char delims[] = "\t, ";
        char *t;
        t = strtok(line, delims);
        edges[offset].from = atoi(t) - 1;
        t = strtok(NULL, delims);
        edges[offset].to = atoi(t) -1 ;

        // if (edges[offset].from == edges[offset].to) continue;

        if ((symmetric == true) && (edges[offset].from != edges[offset].to)) {
            edges[offset + 1].from = edges[offset].to;
            edges[offset + 1].to = edges[offset].from;

            offset++;
        }
        offset++;
    }
    mat.ne = offset;

    sort_To(edges, mat.nr, mat.ne);
    sort_From(edges, mat.nr, mat.ne);
    
    for (int q = 0; q < 10; ++q) {
        free(queues[q]);
    }

    int row_Offset = 0;
    int edges_in_row = 0;
    for (int e = 0; e < mat.ne; ++e) {
        if (edges[e].from > row_Offset) {
            mat.nnz[row_Offset] = edges_in_row;

            for (int empty_Rows = row_Offset + 1; empty_Rows < edges[e].from; ++empty_Rows) {
                mat.nnz[empty_Rows] = 0;
            }

            row_Offset = edges[e].from;
            edges_in_row = 0;
        }
        mat.idx[e] = edges[e].to;
        edges_in_row++;
    }
    mat.nnz[row_Offset] = edges_in_row;

    mat.ptr[0] = 0;
    mat.ptr[mat.nr] = mat.ne;
    for (int row = 0; row < mat.nr - 1; ++row) {
        mat.ptr[row + 1] = mat.ptr[row] + mat.nnz[row];
    }

    free(edges);
    fclose(fp);

    puts("\n --- Input Matrix Information -------------------------------------");
    printf("| Number of Rows:\t%d\n| Number of Edges:\t%d\n\n", mat.nr - 1, mat.ne);

    char input[64];
    strcpy(input, input_File);
    strcat(input, ".csr");
    fp = fopen(input, "wb");
    DATA_TYPE row = mat.nr, col = mat.nc, edge = mat.ne;
    fwrite(&row, sizeof(DATA_TYPE), 1, fp);
    fwrite(&col, sizeof(DATA_TYPE), 1, fp);
    fwrite(&edge, sizeof(DATA_TYPE), 1, fp);

    fwrite(mat.nnz, sizeof(DATA_TYPE), row, fp);
    fwrite(mat.ptr, sizeof(DATA_TYPE), row + 1, fp);
    fwrite(mat.idx, sizeof(DATA_TYPE), edge, fp);
    fclose(fp);

    // free(mat.nnz);
    // free(mat.ptr);
    // free(mat.idx);
}

void convert_synthetic(const char* name) {
    char input[64];
    strcpy(input, name);
    strcat(input, ".csr2");

    // if (0) {
    if (access(input, F_OK) != -1) {
        FILE *fp = fopen(input, "rb");
        fread(&mat.nr, sizeof(DATA_TYPE), 1, fp);
        fread(&mat.nc, sizeof(DATA_TYPE), 1, fp);
        fread(&mat.ne, sizeof(DATA_TYPE), 1, fp);

        mat.nnz = (DATA_TYPE*)malloc(mat.nr * sizeof(DATA_TYPE));
        mat.ptr = (DATA_TYPE*)malloc((mat.nr + 1) * sizeof(DATA_TYPE));
        mat.idx = (DATA_TYPE*)malloc(mat.ne * sizeof(DATA_TYPE));
        fread(mat.nnz, sizeof(DATA_TYPE), mat.nr, fp);
        fread(mat.ptr, sizeof(DATA_TYPE), mat.nr + 1, fp);
        fread(mat.idx, sizeof(DATA_TYPE), mat.ne, fp);

        fclose(fp);

        puts("\n --- Input Matrix Information -------------------------------------");
        printf("| Number of Rows:\t%d\n| Number of Edges:\t%d\n\n", mat.nr - 1, mat.ne);
    }
    else {
        make_synthetic(name);
    }
}


void make_synthetic(const char *input_File) {
    FILE *fp = fopen(input_File, "rt");

    char line[200];
    char *t;
    char delims[] = "\t, ";
    bool symmetric = false;

    fgets(line, 200, fp);
    fgets(line, 200, fp);
    t = strtok(line, delims);
    mat.nr = atoi(t);
    t = strtok(NULL, delims);
    mat.nc = atoi(t);
    t = strtok(NULL, delims);
    mat.ne = atoi(t);


    mat.nnz = (DATA_TYPE*)malloc((mat.nr + 1) * sizeof(DATA_TYPE));
    mat.ptr = (DATA_TYPE*)malloc((mat.nr + 1) * sizeof(DATA_TYPE));
    mat.idx = (DATA_TYPE*)malloc(mat.ne * sizeof(DATA_TYPE));

    edgeData *edges = (edgeData*)malloc(sizeof(edgeData) * mat.ne);
    for (int q = 0; q < 10; ++q) queues[q] = (edgeData*)malloc(sizeof(edgeData) * mat.ne);
    int offset = 0;

    while (fgets(line, 200, fp) != NULL) {
        char delims[] = "\t, ";
        char *t;
        t = strtok(line, delims);
        edges[offset].from = atoi(t) - 1;
        t = strtok(NULL, delims);
        edges[offset].to = atoi(t) -1 ;

        // if (edges[offset].from == edges[offset].to) continue;

        if ((symmetric == true) && (edges[offset].from != edges[offset].to)) {
            edges[offset + 1].from = edges[offset].to;
            edges[offset + 1].to = edges[offset].from;

            offset++;
        }
        offset++;
    }
    mat.ne = offset;

    sort_To(edges, mat.nr, mat.ne);
    sort_From(edges, mat.nr, mat.ne);

    for (int q = 0; q < 10; ++q) {
        free(queues[q]);
    }

    int row_Offset = 0;
    int edges_in_row = 0;
    offset = 0;
    for (int e = 0; e < mat.ne; ++e) {
        if (edges[e].from > row_Offset) {
            mat.nnz[row_Offset] = edges_in_row;

            for (int empty_Rows = row_Offset + 1; empty_Rows < edges[e].from; ++empty_Rows) {
                mat.nnz[empty_Rows] = 0;
            }

            row_Offset = edges[e].from;
            edges_in_row = 0;
        }
        if (edges[e].to == edges[e - 1].to) continue;
        // if (edges[e].to <= edges[e - 1].to) edges[e].to = edges[e - 1].to + 1;
        // if (edges[e].to >= mat.nr) continue;
        mat.idx[e] = edges[e].to;
        edges_in_row++;
    }
    mat.nnz[row_Offset] = edges_in_row;

    mat.ptr[0] = 0;
    // mat.ptr[mat.nr] = mat.ne;
    for (int row = 0; row < mat.nr; ++row) {
        mat.ptr[row + 1] = mat.ptr[row] + mat.nnz[row];
    }
    mat.ne = mat.ptr[mat.nr];

    free(edges);
    fclose(fp);

    puts("\n --- Input Matrix Information -------------------------------------");
    printf("| Number of Rows:\t%d\n| Number of Edges:\t%d\n\n", mat.nr - 1, mat.ne);

    char input[64];
    strcpy(input, input_File);
    strcat(input, ".csr2");
    fp = fopen(input, "wb");
    DATA_TYPE row = mat.nr, col = mat.nc, edge = mat.ne;
    fwrite(&row, sizeof(DATA_TYPE), 1, fp);
    fwrite(&col, sizeof(DATA_TYPE), 1, fp);
    fwrite(&edge, sizeof(DATA_TYPE), 1, fp);

    fwrite(mat.nnz, sizeof(DATA_TYPE), row, fp);
    fwrite(mat.ptr, sizeof(DATA_TYPE), row + 1, fp);
    fwrite(mat.idx, sizeof(DATA_TYPE), edge, fp);
    fclose(fp);

    // free(mat.nnz);
    // free(mat.ptr);
    // free(mat.idx);
}


/** 
 * Radix Sort
 * Each edge has {from, to, value} data
 * Sort only with from data
 * Total Memory Usage = sizeof(struct edgeData) * mat.ne +                                               : Temporary Array(edges in main Function) to read input file
 *                      10 * sizeof(struct edgeData) * mat.ne +                                                      : Queue Size
 *                      sizeof(int) * mat.nr + sizeof(int) * mat.ne + sizeof(float) * mat.ne        : Output Array(CSR Format)
 *                      Actually, Ouput Array Size is hidden since memory allocated after freeing Queues
 * Can reduce memory usage by changing implemention to 2-bit radix sort
 * -> Queue Size = 2 * sizeof(int) * mat.ne
 * But, time complexity increases.
 */

void sort_From(struct edgeData *edges, int nr, int ne) {
    int num_Digits = 1;
    while (nr) {
        nr /= 10;
        num_Digits *= 10;
    }
    //printf("Number of Digits: %d\n", num_Digits);

    for (int divide_To = 1; divide_To <= num_Digits; divide_To *= 10) {
        for (int q = 0; q < 10; ++q) queue_Offsets[q] = 0;
        
        for (int e = 0; e < ne; ++e) {
            int row = edges[e].from;
            int target_Queue = row % divide_To;
            target_Queue = (target_Queue * 10) / divide_To;
            queues[target_Queue][queue_Offsets[target_Queue]++] = edges[e];
        }

        int offset = 0;
        for (int q = 0; q < 10; ++q) {
            memcpy(edges + offset, queues[q], sizeof(struct edgeData) * queue_Offsets[q]);
            offset += queue_Offsets[q];
        }
    }
}

void sort_To(struct edgeData *edges, int nc, int ne) {
    int num_Digits = 1;
    while (nc) {
        nc /= 10;
        num_Digits *= 10;
    }
    //printf("Number of Digits: %d\n", num_Digits);

    for (int divide_To = 1; divide_To <= num_Digits; divide_To *= 10) {
        for (int q = 0; q < 10; ++q) queue_Offsets[q] = 0;
        
        for (int e = 0; e < ne; ++e) {
            int col = edges[e].to;
            int target_Queue = col % divide_To;
            target_Queue = (target_Queue * 10) / divide_To;
            queues[target_Queue][queue_Offsets[target_Queue]++] = edges[e];
        }

        int offset = 0;
        for (int q = 0; q < 10; ++q) {
            memcpy(edges + offset, queues[q], sizeof(struct edgeData) * queue_Offsets[q]);
            offset += queue_Offsets[q];
        }
    }
}
