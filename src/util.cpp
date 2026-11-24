#include "SpMM.hpp"

struct time_log {
    char str[64];
    float time = 0.0f;
    timeval st, ed;
};

typedef struct time_log timelog;
typedef struct timeval timeval;

#define NUM_LOG 9

timelog tl[NUM_LOG];
int fd_dma;
DATA_TYPE total_merge = 0, total_split = 0, total_segment = 0, total_gemm = 0, total_binded_tile = 0;
float total_time = 0.0f;

inline float elapsed_Time(timeval st, timeval ed) {return 1000.0f * (ed.tv_sec - st.tv_sec) + 0.001f * (ed.tv_usec - st.tv_usec);}
 
void cur_Time(timeval *tv) {
    gettimeofday(tv, NULL);
}

void start_Time(enum TIME_LOG type) {
    cur_Time(&tl[type].st);
}

void end_Time(enum TIME_LOG type) {
    cur_Time(&tl[type].ed);
    tl[type].time += elapsed_Time(tl[type].st, tl[type].ed);
}

void print_TimeLog(int mode) {
    strcpy(tl[TIME_INIT].str, "Initialize Input Matrices");
    strcpy(tl[TIME_TCSR].str, "Build Tiled CSR Format");
    strcpy(tl[TIME_MERGE].str, "Merge Identical Segments");
    strcpy(tl[TIME_ALIGN].str, "Align Input Matrices for PIM");
    strcpy(tl[TIME_DMA].str, "DMA Tranfer from Host to PIM");
    strcpy(tl[TIME_MATMUL].str, "Perform SpMM with GEMM APIs");
    strcpy(tl[TIME_OVERLAP].str, "Overlapped SpMM + PIM Align");
    strcpy(tl[TIME_TEST].str, "TEST 1");
    strcpy(tl[TIME_TEST2].str, "TEST 2");


    puts(" --- Time ---------------------------------------------------------");
    for (int i = 0; i < NUM_LOG; ++i) {
        if (tl[i].time > 0.001f) {
            tl[i].time /= ITER;
            printf("| %s: %.3f msec\n", tl[i].str, tl[i].time);
            if (i < 7) total_time += tl[i].time;
        } 
    }
    if (check_option(mode, OVERLAPPING)) total_time = total_time - tl[TIME_MERGE].time - tl[TIME_ALIGN].time - tl[TIME_DMA].time - tl[TIME_MATMUL].time;
    printf("| Total Time: %.3f msec\n\n", total_time);
    printf("GFLOPS:\n%lf\n", 2.0 * mat.ne * COL_B / ((double)total_time * 1e6));
}

float print_Result(int mode) {
    if (check_option(mode, TCSR_SPMM)) {
        puts(" --- Size of CSR / tCSR Format ------------------------------------");
        printf("| Total Size of CSR Format: %u KB\n", (mat.nr + mat.ne) * sizeof(DATA_TYPE) >> 10);
        printf("| Total Size of tCSR Format: %u KB (%u Segments)\n\n", tmat.nr * sizeof(DATA_TYPE) + tmat.ne * (sizeof(uint16_t) + sizeof(DATA_TYPE)) >> 10, tmat.ne);
    }

    // printf("Total Size of CSR Format: %d KB\n", (mat.nr + mat.ne) * sizeof(int) >> 10);
    // printf("Total Size of tCSR Format: %d KB (%d Segments)\n", (tmat.nr * sizeof(int) + (sizeof(int) + sizeof(uint16_t)) * tmat.ne) >> 10, tmat.ne);
    print_TimeLog(mode);

    if (check_option(mode, DEBUGGING)) {
        int csr_size = (mat.nr + mat.ne) * sizeof(DATA_TYPE) >> 10;
        int tcsr_size = 0;

        FILE *fp = fopen("SDGEMM_RESULTS.csv", "at");
        for (int i = 1; i < 7; ++i) {
            if (tl[i].time > 0.001f) fprintf(fp, ",%f", tl[i].time);
            else fprintf(fp, ",0");
        }
        if (check_option(mode, TCSR_SPMM)) {
            if (!check_option(mode, MERGING)) {
                total_segment = tmat.ne;
                for (int i = 0; i < tmat.nr; ++i) {
                    total_gemm += (tmat.ptr[i + 1] - tmat.ptr[i] + 15) >> 4;
                }
            }
            tcsr_size = tmat.nr * sizeof(DATA_TYPE) + tmat.ne * (sizeof(uint16_t) + sizeof(DATA_TYPE)) >> 10;
        }
        fprintf(fp, ",%f,%d,%d,%d,%d,%d,%d,%d,%d\n", total_time, total_merge, total_split, total_segment, total_gemm, tmat.nr, total_binded_tile, csr_size, tcsr_size);
        fclose(fp);
    }
    return total_time;
}

void print_Merge(workload *wl) {
    int total_merge = 0, total_split = 0, total_segment = 0, total_gemm = 0;
    for (int tile = 0; tile < tmat.nr; ++tile) {
        total_merge += wl[tile].num_Merge;
        total_split += wl[tile].num_Split;
        total_segment += wl[tile].num_Segment;
        total_gemm += wl[tile].num_GEMM;
        // printf("%d %d %d %d\n", total_merge, total_split, total_segment, total_gemm);
    }
    puts(" --- Merge Identical Segments -------------------------------------");
    printf("| Number of Merge: %d\n", total_merge);
    printf("| Number of Split: %d\n", total_split);
    printf("| Number of Segments: %d\n", total_segment);
    printf("| Number of GEMM Units: %d\n\n", total_gemm);
}

void print_Binding(int binded_tile) {
    total_binded_tile = binded_tile;
    puts(" --- Bind Small Tiles with Fused Kernel ---------------------------");
    printf("| Number of Original Tile: %d\n", tmat.nr);
    printf("| Number of Binded Tile: %d\n\n", binded_tile);
}


// void transfer_PIM(short *src, short *dst, DATA_TYPE size) {
//     if (size == 0) return;
//     start_Time(TIME_DMA);
//     pim_args xfer;
//     xfer.srcA_ptr = src;
//     xfer.srcB_ptr = NULL;
//     xfer.dstC_ptr = dst;
//     xfer.srcA_va = (uint64_t)src;
//     xfer.srcB_va = 0x0;
//     xfer.dstC_va = (uint64_t)dst;
//     xfer.srcA_size = size;
//     xfer.srcB_size = 0x0;
//     xfer.dstC_size = size;
//     xfer.p_size = 0x0;
//     xfer.q_size = 0x0;
//     xfer.r_size = 0x0;

//     if (ioctl(fd_dma, MEMCPY_PS2PL, &xfer) < 0) {
//         printf("ERROR DMA \n");
//         return;
//     }
//     end_Time(TIME_DMA);
// }

void matmul_cpu(VAL_TYPE *A, VAL_TYPE *B, VAL_TYPE *C, int padded_segment) {
    start_Time(TIME_MATMUL);
// #pragma omp parallel for num_threads(4) schedule(static, 4) 
#pragma omp parallel for num_threads(4) 
    for (int i = 0; i < NUM_BANK; ++i) {
        for (int k = 0; k < padded_segment; ++k) {
            for (int j = 0; j < COL_B; ++j) {
                C[i * COL_B + j] += A[i * padded_segment + k] * B[k * COL_B + j];
            }
        }
    }
    end_Time(TIME_MATMUL);
}

void binded_matmul_cpu(VAL_TYPE *A, VAL_TYPE *B, VAL_TYPE *C, uint32_t *num_segment, uint32_t num_tile) {
    start_Time(TIME_MATMUL);
    VAL_TYPE *A_off = A, *B_off = B, *C_off = C;
    for (int tile = 0; tile < num_tile; ++tile) {
        for (int i = 0; i < NUM_BANK; ++i) {
            for (int k = 0; k < num_segment[tile]; ++k) {
                for (int j = 0; j < COL_B; ++j) {
                    C_off[i * COL_B + j] += A_off[i * num_segment[tile] + k] * B_off[k * COL_B + j];
                }
            }
        }

        A_off += NUM_BANK * num_segment[tile];
        B_off += num_segment[tile] * COL_B;
        C_off += COL_B * NUM_BANK;
    }   
    end_Time(TIME_MATMUL);
}

int emul_time[15] = {25, 37, 49, 61, 76, 88, 100, 112, 126, 138, 150, 165, 177, 190, 202};
#define SCALED_TIME_PER_16_SEGMENTS 12.911
void emulate_gemm_pim(VAL_TYPE *A, VAL_TYPE *B, VAL_TYPE *C, uint32_t padded_segment) {
    start_Time(TIME_MATMUL);
    int usec;
    if (padded_segment < 256) usec = emul_time[padded_segment >> 4];
    else usec = (int)round((padded_segment >> 4) * SCALED_TIME_PER_16_SEGMENTS);

    usleep(usec);
    end_Time(TIME_MATMUL);
}

int emul_perbank_time[15] = {34, 58, 79, 104, 125, 148, 170, 197, 214, 236, 257, 279, 301, 323, 345};
#define SCALED_TIME_PER_16_SEGMENTS 21.991
void emulate_perbank_pim(VAL_TYPE *A, VAL_TYPE *B, VAL_TYPE *C, uint32_t padded_segment) {
    start_Time(TIME_MATMUL);
    int usec;
    if (padded_segment < 256) usec = emul_perbank_time[padded_segment >> 4];
    else usec = (int)round((padded_segment >> 4) * SCALED_TIME_PER_16_SEGMENTS);

    usleep(usec);
    end_Time(TIME_MATMUL);
}

void analyze_Locality() {
    char *dense = (char*)malloc(mat.nc);
    char *flag = (char*)calloc(mat.ne, 1);

    DATA_TYPE count = 0;
    for (int i = 0; i < mat.nr; ++i) {
        for (int element = mat.ptr[i]; element < mat.ptr[i + 1]; ++element) {
            dense[mat.idx[element]] = 1;
        }
        if (i > 0) {
            for (int element = mat.ptr[i - 1]; element < mat.ptr[i]; ++element) {
                if (dense[mat.idx[element]] == 1) flag[element] = 1;
            }
        }
        if (i < mat.nr - 1) {
            for (int element = mat.ptr[i + 1]; element < mat.ptr[i + 2]; ++element) {
                if (dense[mat.idx[element]] == 1) flag[element] = 1;
            }
        }

        for (int element = mat.ptr[i]; element < mat.ptr[i + 1]; ++element) {
            dense[mat.idx[element]] = 0;
        }
    }

    for (int i = 0; i < mat.ne; ++i) {
        if (flag[i]) count++;
    }

    printf("Ratio of non-zeros with data locality: %f\n", (float)count / mat.ne);
}

void sort_workload(workload* wl, int* sorted, int mode) {
    DATA_TYPE *arr = (DATA_TYPE*)malloc(tmat.nr << SHIFT_BIT);
    if (check_option(mode, MERGING)) {
        for (int i = 0; i < tmat.nr; ++i) {
            arr[i] = wl[i].num_GEMM;
        }
    }
    else {
        for (int i = 0; i < tmat.nr; ++i) {
            arr[i] = tmat.ptr[i + 1] - tmat.ptr[i];
        }
    }

    for (int tile = 0; tile < tmat.nr; ++tile) {
        int max = 0xffffffff;
        int idx;
        for (int i = 0; i < tmat.nr; ++i) {
            if (arr[i] > max) {
                max = arr[i];
                idx = i;
            }
        }
        arr[idx] = 0;
        sorted[tile] = idx;
    }
}
