#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/ioctl.h>
#include <time.h>
#include <sys/time.h>
#include <assert.h>
#include <limits.h>
#include <iostream>
#include <math.h>
#include <bitset>
#include <pthread.h>
#include <omp.h>

#define NUM_BANK 16
#define REG_SIZE 16
#define COL_B 1024
enum TIME_LOG {TIME_INIT, TIME_TCSR, TIME_MERGE, TIME_ALIGN, TIME_DMA, TIME_MATMUL, TIME_OVERLAP, TIME_TEST, TIME_TEST2};
#define BASE_SPMM 1 // CPU-implemented SpMM
#define UNCOMPRESSED 2
#define CSR_SPMM 4
#define TCSR_SPMM 8
#define MERGING 16
#define OVERLAPPING 32
#define BINDING 64
#define PARALLEL 128
#define DEBUGGING 256

#define DATA_TYPE int32_t
#define SHIFT_BIT 2
#define VAL_TYPE short

#define ITER 3
// #define DATA_TYPE uint16_t
// #define SHIFT_BIT 1

#define check_option(mode, option) ((mode & option) != 0)
#define PADDING(num) ((num + 15) >> 4) << 4

#define THRESHOLD 8

struct CSR {
    DATA_TYPE *nnz;
    DATA_TYPE *ptr;
    DATA_TYPE *idx;
    DATA_TYPE *mapping; // Row mapping
    DATA_TYPE nr;
    DATA_TYPE nc;
    DATA_TYPE ne;
};
typedef struct CSR csr;


struct tCSR {
    DATA_TYPE* ptr;
    // DATA_TYPE* idx; 
    DATA_TYPE* mapping; // Segment mapping
    // DATA_TYPE* col_from;
    uint8_t* num;
    uint16_t* row_hash;
    DATA_TYPE nr;
    DATA_TYPE ne;
};
typedef struct tCSR tcsr;

struct MERGE {
    DATA_TYPE src; // Actual row
    DATA_TYPE dst; // Segment index
};
typedef struct MERGE merge;

struct SPLIT {
    DATA_TYPE src;
    DATA_TYPE dst1;
    DATA_TYPE dst2;
};
typedef struct SPLIT split;

struct WORKLOAD {
    DATA_TYPE num_Segment;
    DATA_TYPE num_Merge;
    DATA_TYPE num_GEMM;
    DATA_TYPE num_Split;
};
typedef struct WORKLOAD workload;

// void analyze_Locality();

// Input Matrix Preprocessing
void convert_CSR(const char *name);
void convert_synthetic(const char *name);
void make_TiledCSR(int mode);

// For Segment Merge
void init_Merge(int mode);
void end_Merge(int mode);
void merge_Segments(workload *wl, merge *merging, split *splitting);
void merge_Tile(workload *wl, merge *merging, split *splitting, int tidx);
void merge_Tile_Non_Binary(workload *wl, merge *merging, split *splitting, int tidx);

// For Performance Comparison
void SpMM_CPU();
void GEMM_CPU();
void tiled_spmm(VAL_TYPE *align_A, VAL_TYPE *B, VAL_TYPE *align_B, VAL_TYPE *C, uint32_t M, uint32_t K, uint32_t N, int mode);

// For Matrix Multiplication
void matmul_cpu(VAL_TYPE *A, VAL_TYPE *B, VAL_TYPE *C, int padded_segment);
void binded_matmul_cpu(VAL_TYPE *A, VAL_TYPE *B, VAL_TYPE *C, uint32_t *num_segment, uint32_t num_tile);
void emulate_gemm_pim(VAL_TYPE *A, VAL_TYPE *B, VAL_TYPE *C, uint32_t padded_segment);
void emulate_perbank_pim(VAL_TYPE *A, VAL_TYPE *B, VAL_TYPE *C, uint32_t padded_segment);

// For Data Transfer with DMA
void transfer_PIM(VAL_TYPE *src, VAL_TYPE *dst, DATA_TYPE size);

// Supporting Functions for PIM
void rand_mat(short *m, int R, int S);
short float_to_short(float value);
float short_to_float (short value);

// For Printing and Analyzing Results
float print_Result(int mode);
void print_Merge(workload *wl);
void print_Binding(int binded_tile);

void start_Time(enum TIME_LOG type);
void end_Time(enum TIME_LOG type);


extern csr mat;
extern tcsr tmat;
extern int fd_dma;
