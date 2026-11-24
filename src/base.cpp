#include "SpMM.hpp"
 
// float *val;
// float *B;
// float *C;

// void *spmm_parallel(void *args) {
//     for (int row = 0; row < mat.nr; ++row) {
//         for (int element = mat.ptr[row]; element < mat.ptr[row + 1]; ++element) {
//             int row_b = mat.idx[element];
//             for (int col_b = (col_b >> 1); col_b < COL_B; ++col_b) {
//                 C[row * COL_B + col_b] += val[element] * B[row_b * COL_B + col_b];
//             }
//         }
//     }

//     return NULL;
// }
 

// Row Product
void SpMM_CPU() {
    struct timeval st, ed;
    DATA_TYPE M = mat.nr;
    DATA_TYPE K = mat.nc;
    DATA_TYPE N = COL_B;
    float *val = (float*)malloc(mat.ne * sizeof(float));
    float *B = (float*)malloc(K * N * sizeof(float));
    float *C = (float*)malloc(M * N * sizeof(float));
    for (int i = 0; i < mat.ne; ++i) val[i] = 1.0f;
    for (int i = 0; i < K * N; ++i) B[i] = 3.0f;
    memset(C, 0, M * N << 2);
    
    start_Time(TIME_MATMUL);
    gettimeofday(&st, NULL);
    // pthread_t worker;
    // pthread_create(&worker, NULL, SpMM_parallel, NULL);

    for (int row_a = 0; row_a < mat.nr; ++row_a) {
        for (int element = mat.ptr[row_a]; element < mat.ptr[row_a + 1]; ++element) {
            int col_a = mat.idx[element];
            for (int col_b = 0; col_b < N; ++col_b) {
                C[row_a * N + col_b] += val[element] * B[col_a * N + col_b];
            }
        }
    }
    // pthread_join(worker, NULL);
    gettimeofday(&ed, NULL);
    double time = 1e6 * (ed.tv_sec - st.tv_sec) + 1.0 * (ed.tv_usec - st.tv_usec);
    printf("\nSpMM:\n%lf\n", time * 1e-3);
    printf("%lf\n", 2.0 * mat.ne * COL_B / time / 1e3) ;
    end_Time(TIME_MATMUL);

    float acc = 0.0f;
    for (int i = 0; i < M * N; ++i) acc += C[i];
    // printf("Acc: %f\n", acc);

    free(B);
    free(C);
    free(val);
}


// void *spmm_parallel2(void *args) {
//     for (int row_a = 0; row_a < mat.nr; ++row_a) {
//         for (int col_b = 0; col_b < (COL_B >> 1); ++col_b) {
//             for (int element = mat.ptr[row_a]; element < mat.ptr[row_a + 1]; ++element) {
//                 int col_a = mat.idx[element];
//                 C[row_a * COL_B + col_b] += val[element] * B[col_a * COL_B + col_b];
//             }
//         }
//     }

//     return NULL;
// }

// Inner Product
void SpMM_CPU2() {
    struct timeval st, ed;
    DATA_TYPE M = mat.nr;
    DATA_TYPE K = mat.nc;
    DATA_TYPE N = COL_B;
    float *val = (float*)malloc(mat.ne << 2);
    for (int i = 0; i < mat.ne; ++i) {
        val[i] = 1.0f;
    }

    float *B = (float*)malloc(K * N << 2);
    for (int i = 0; i < K * N; ++i) {
        B[i] = 3.0f;
    }
    float *C = (float*)malloc(M * N << 2);
    memset(C, 0, M * N << 2);
    
    start_Time(TIME_MATMUL);
    gettimeofday(&st, NULL);
    // pthread_t worker;
    // pthread_create(&worker, NULL, spmm_parallel2, NULL);

    for (int row_a = 0; row_a < mat.nr; ++row_a) {
        for (int col_b = 0; col_b < (COL_B >> 1); ++col_b) {
            for (int element = mat.ptr[row_a]; element < mat.ptr[row_a + 1]; ++element) {
                int col_a = mat.idx[element];
                C[row_a * N + col_b] += val[element] * B[col_a * N + col_b];
            }
        }
    }
    // pthread_join(worker, NULL);
    gettimeofday(&ed, NULL);
    double time = 1000000.0 * (ed.tv_sec - st.tv_sec) + (double)(ed.tv_usec - st.tv_usec);
    printf("Inner Product: %.3lf\n", time / 1000.0);
    printf("GFLOPS: %lf\n", (double)2 * mat.ne * COL_B / time / 1000.0) ;
    end_Time(TIME_MATMUL);
    free(B);
    free(C);
    free(val);
}

// GFLOPS Test
void GEMM_CPU() {
    struct timeval st, ed;
    int size = (int)ceil(sqrt(mat.ne));
    int M = size;
    int K = size;
    int N = COL_B;
    float *A = (float*)malloc(size * size << 2);
    float *B = (float*)malloc(size * COL_B << 2);
    float *C = (float*)malloc(size * COL_B << 2);
    for (int i = 0; i < M * K; ++i) A[i] = 1.0f;
    for (int i = 0; i < K * N; ++i) B[i] = 3.0f;
    memset(C, 0, M * N << 2);

    
    gettimeofday(&st, NULL);
    // matmul_cpu(A, B, C, M, N, K);
    for (int row_a = 0; row_a < M; ++row_a) {
        // printf("Row %d\n", row_a);
        for (int col_a = 0; col_a < K; ++col_a) {
            for (int col_b = 0; col_b < N; ++col_b) {
                C[row_a * N + col_b] += A[row_a * K + col_a] * B[col_a * N + col_b];
            }

        }
    }
    gettimeofday(&ed, NULL);
    // double time = 1000000.0 * (ed.tv_sec - st.tv_sec) + (double)(ed.tv_usec - st.tv_usec);
    double time = 1e6 * (ed.tv_sec - st.tv_sec) + 1.0 * (ed.tv_usec - st.tv_usec);
    printf("\nGEMM:\n%lf\n", time * 1e-3);
    printf("%lf\n", 2.0 * size * size * COL_B / time / 1e3);

    float acc = 0.0f;
    for (int i = 0; i < M * N; ++i) acc += C[i];
    printf("Acc: %f\n", acc);

    free(A);
    free(B);
    free(C);
}

// void GEMM_CPU2() {
//     struct timeval st, ed;
//     int size = (int)ceil(sqrt((float)mat.ne));
//     float *A = (float*)malloc(size * size << 2);
//     float *B = (float*)malloc(size * COL_B << 2);
//     float *C = (float*)malloc(size * COL_B << 2);

//     for (int i = 0; i < size * size; ++i) A[i] = 1.0f;
//     for (int i = 0; i < size * COL_B; ++i) A[i] = 3.0f;
//     memset(C, 0, size * COL_B << 2);

    
//     gettimeofday(&st, NULL);
//     for (int row_a = 0; row_a < size; ++row_a) {
//         for (int col_b = 0; col_b < COL_B; ++col_b) {
//             for (int col_a = 0; col_a < size; ++col_a) {
//                 C[row_a * COL_B + col_b] += A[row_a * size + col_a] * B[col_a * COL_B + col_b];
//             }
//         }
//     }
//     gettimeofday(&ed, NULL);
//     double time = 1000000.0 * (ed.tv_sec - st.tv_sec) + (double)(ed.tv_usec - st.tv_usec);
//     printf("Time: %lf\n", time * 0.001);
//     printf("%lf\n", (double)2 * size * size * COL_B / time / 1000.0) ;

//     free(A);
//     free(B);
//     free(C);
// }