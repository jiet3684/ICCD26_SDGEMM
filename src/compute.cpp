#include "SpMM.hpp"

#define ONE 16256
#define BIND_CRITERIA 512

int status;
#define WAITING 0
#define WORKING 1

merge *merging;
split *splitting;
VAL_TYPE *merge_split;
workload *wl;
uint8_t *moved;

DATA_TYPE *col_to;
DATA_TYPE *col_from;

void init_SpMM(int mode) {
    if (check_option(mode, MERGING)) {
        init_Merge(mode);
        merging = (merge*)malloc(tmat.ne * sizeof(merge));
        splitting = (split*)malloc(tmat.nr * THRESHOLD * sizeof(split));
        merge_split = (VAL_TYPE*)malloc(THRESHOLD * COL_B * NUM_BANK * sizeof(VAL_TYPE));
        memset(merge_split, 0, THRESHOLD * COL_B * sizeof(VAL_TYPE));
        wl = (workload*)malloc(tmat.nr * sizeof(workload));
        // merge_Segments(wl, merging, splitting);
    }
}

void end_SpMM(int mode) {
    if (check_option(mode, MERGING)) {
        print_Merge(wl);

        free(merging);
        free(splitting);
        free(merge_split);
        free(wl);
        end_Merge(mode);
    }
}
uint32_t align_tile_parallel(VAL_TYPE *align_A, VAL_TYPE *B, VAL_TYPE *align_B, int tidx, int mode) {
    start_Time(TIME_ALIGN);

    int num_segment = tmat.ptr[tidx + 1] - tmat.ptr[tidx];
    if (check_option(mode, MERGING)) num_segment = wl[tidx].num_Segment;
    int num_gemm = (num_segment + NUM_BANK - 1) / NUM_BANK; 
    int padded_segment = num_gemm * NUM_BANK;
    uint32_t tile_offset = tmat.ptr[tidx];
    int cnt = 0;

// Sparse matrix A align
// #pragma omp parallel for num_threads(4) 
    for (uint32_t segment = tile_offset; segment < tile_offset + num_segment; ++segment) {
        uint16_t hash_val = tmat.row_hash[segment];
        int offset = (segment - tile_offset) << 4;
        align_A[offset    ] = hash_val & 0x0001;
        align_A[offset + 1] = hash_val & 0x0002 >> 1;
        align_A[offset + 2] = hash_val & 0x0004 >> 2;
        align_A[offset + 3] = hash_val & 0x0008 >> 3;
        align_A[offset + 4] = hash_val & 0x0010 >> 4;
        align_A[offset + 5] = hash_val & 0x0020 >> 5;
        align_A[offset + 6] = hash_val & 0x0040 >> 6;
        align_A[offset + 7] = hash_val & 0x0080 >> 7;
        align_A[offset + 8] = hash_val & 0x0100 >> 8;
        align_A[offset + 9] = hash_val & 0x0200 >> 9;
        align_A[offset + 10] = hash_val & 0x0400 >> 10;
        align_A[offset + 11] = hash_val & 0x0800 >> 11;
        align_A[offset + 12] = hash_val & 0x1000 >> 12;
        align_A[offset + 13] = hash_val & 0x2000 >> 13;
        align_A[offset + 14] = hash_val & 0x4000 >> 14;
        align_A[offset + 15] = hash_val & 0x8000 >> 15;
    }
    // for (uint32_t segment = tile_offset; segment < tile_offset + num_segment; ++segment) {
    //     uint16_t hash_val = tmat.row_hash[segment];
    //     align_A[(cnt << 4)    ] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 1] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 2] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 3] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 4] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 5] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 6] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 7] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 8] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 9] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 10] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 11] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 12] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 13] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 14] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     align_A[(cnt << 4) + 15] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
    //     cnt++;
    // }
    
    for (uint32_t pad = num_segment; pad < padded_segment; ++pad) {
        for (int element = 0; element < 16; ++element) align_A[element * padded_segment + pad] = 0;
    }

// Weight matrix B align
    cnt = 0;
// #pragma omp parallel for num_threads(4)
    for (int row = 0; row < num_segment; ++row) {
        memcpy(align_B + cnt, B + (tmat.mapping[tile_offset + row] * COL_B), COL_B * sizeof(VAL_TYPE));
        cnt += COL_B;
    }
    memset(align_B + cnt, 0, (padded_segment - num_segment) * COL_B * sizeof(VAL_TYPE));
    #define NUM_THREAD 4


    start_Time(TIME_TEST);
    if (check_option(mode, MERGING)) {
        memset(merge_split, 0x00, THRESHOLD * COL_B * NUM_BANK * sizeof(VAL_TYPE));
// Merge corresponding rows of matrix B
#pragma omp parallel for num_threads(NUM_THREAD)
        for (int segment = tile_offset; segment < tile_offset + wl[tidx].num_Merge; ++segment) {
            VAL_TYPE *my_merge_split = merge_split + (omp_get_thread_num() * THRESHOLD * COL_B);
            DATA_TYPE source_row = merging[segment].src;
            DATA_TYPE target_segment = merging[segment].dst;
            if ((target_segment & 0x80000000) == 0) {
                DATA_TYPE dst = target_segment * COL_B, src = source_row * COL_B;
#pragma unroll(8)
                for (int i = 0; i < COL_B; ++i) {
                    // align_B[target_segment * COL_B + i] += B[(source_row * COL_B) + i];
                    align_B[dst + i] += B[src + i];
                }
            }
            else {
                DATA_TYPE target = merging[segment].dst & 0x0000000f;
                for (int i = 0; i < COL_B; ++i) {
                    my_merge_split[(target * COL_B) + i] += B[(source_row * COL_B) + i];
                }
            }
            // if ((target_segment & 0x80000000) == 0x80000000) {
            //     DATA_TYPE target = merging[segment].dst & 0x0000000f;
            //     for (int i = 0; i < COL_B; ++i) {
            //         merge_split[(target * COL_B) + i] += B[(source_row * COL_B) + i];
            //     }
            // }
            // else {
            //     for (int i = 0; i < COL_B; ++i) {
            //         align_B[target_segment * COL_B + i] += B[(source_row * COL_B) + i];
            //     }
            // }
        }
        for (int thread = 1; thread < NUM_THREAD; ++thread) {
            for (int col = 0; col < COL_B; ++col) {
                merge_split[col] += merge_split[thread * COL_B + col];
            }
        }
        end_Time(TIME_TEST);

        start_Time(TIME_TEST2);
// Split corresponding rows of matrix B
        for (int segment = tidx * THRESHOLD; segment < tidx * THRESHOLD + wl[tidx].num_Split; ++segment) {
            DATA_TYPE source_row = splitting[segment].src;
            DATA_TYPE target_segment1 = splitting[segment].dst1;
            DATA_TYPE target_segment2 = splitting[segment].dst2;
            VAL_TYPE *source;
            if (source_row & 0x80000000) {
                DATA_TYPE target = splitting[segment].src & 0x0000000f;
                source = merge_split + (target * COL_B);
            }
            else source = B + (source_row * COL_B);
            for (int i = 0; i < COL_B; ++i) {
                align_B[target_segment1 * COL_B + i] += source[i];
                align_B[target_segment2 * COL_B + i] += source[i];
            }
        }

        // for (int segment = tile_offset; segment < tile_offset + wl[tidx].num_Merge; ++segment) {
        //     if (merging[segment].dst & 0x80000000) {
        //         DATA_TYPE target = merging[segment].dst & 0x0000000f;
        //         memset(merge_split + (target * COL_B), 0, COL_B * sizeof(float));
        //     }
        // }
    }
    end_Time(TIME_TEST2);

    end_Time(TIME_ALIGN);

    return padded_segment;
}


uint32_t align_tile(VAL_TYPE *align_A, VAL_TYPE *B, VAL_TYPE *align_B, int tidx, int mode) {
    start_Time(TIME_ALIGN);

    int num_segment = tmat.ptr[tidx + 1] - tmat.ptr[tidx];
    if (check_option(mode, MERGING)) num_segment = wl[tidx].num_Segment;
    int num_gemm = (num_segment + NUM_BANK - 1) / NUM_BANK; 
    int padded_segment = num_gemm * NUM_BANK;
    uint32_t tile_offset = tmat.ptr[tidx];
    int cnt = 0;

// Sparse matrix A align
    for (uint32_t segment = tile_offset; segment < tile_offset + num_segment; ++segment) {
        uint16_t hash_val = tmat.row_hash[segment];
        align_A[(cnt << 4)    ] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 1] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 2] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 3] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 4] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 5] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 6] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 7] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 8] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 9] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 10] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 11] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 12] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 13] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 14] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        align_A[(cnt << 4) + 15] = hash_val & 0x0001 == 0x0000 ? 0.0f : 1.0f;
        cnt++;
    }
    
    for (uint32_t pad = num_segment; pad < padded_segment; ++pad) {
        for (int element = 0; element < 16; ++element) align_A[element * padded_segment + pad] = 0;
    }

// Weight matrix B align
    cnt = 0;
    for (int row = 0; row < num_segment; ++row) {
        memcpy(align_B + cnt, B + (tmat.mapping[tile_offset + row] * COL_B), COL_B * sizeof(VAL_TYPE));
        cnt += COL_B;
    }
    memset(align_B + cnt, 0, (padded_segment - num_segment) * COL_B * sizeof(VAL_TYPE));


    if (check_option(mode, MERGING)) {
        start_Time(TIME_TEST);
        memset(merge_split, 0x00, THRESHOLD * COL_B * sizeof(VAL_TYPE));
// Merge corresponding rows of matrix B
        for (int segment = tile_offset; segment < tile_offset + wl[tidx].num_Merge; ++segment) {
            DATA_TYPE source_row = merging[segment].src;
            DATA_TYPE target_segment = merging[segment].dst;
            if ((target_segment & 0x80000000) == 0) {
                for (int i = 0; i < COL_B; ++i) {
                    align_B[target_segment * COL_B + i] += B[(source_row * COL_B) + i];
                }
            }
            else {
                DATA_TYPE target = merging[segment].dst & 0x0000000f;
                for (int i = 0; i < COL_B; ++i) {
                    merge_split[(target * COL_B) + i] += B[(source_row * COL_B) + i];
                }
            }
            // if ((target_segment & 0x80000000) == 0x80000000) {
            //     DATA_TYPE target = merging[segment].dst & 0x0000000f;
            //     for (int i = 0; i < COL_B; ++i) {
            //         merge_split[(target * COL_B) + i] += B[(source_row * COL_B) + i];
            //     }
            // }
            // else {
            //     for (int i = 0; i < COL_B; ++i) {
            //         align_B[target_segment * COL_B + i] += B[(source_row * COL_B) + i];
            //     }
            // }
        }
        end_Time(TIME_TEST);

        start_Time(TIME_TEST2);
// Split corresponding rows of matrix B
        for (int segment = tidx * THRESHOLD; segment < tidx * THRESHOLD + wl[tidx].num_Split; ++segment) {
            DATA_TYPE source_row = splitting[segment].src;
            DATA_TYPE target_segment1 = splitting[segment].dst1;
            DATA_TYPE target_segment2 = splitting[segment].dst2;
            VAL_TYPE *source;
            if (source_row & 0x80000000) {
                DATA_TYPE target = splitting[segment].src & 0x0000000f;
                source = merge_split + (target * COL_B);
            }
            else source = B + (source_row * COL_B);
            for (int i = 0; i < COL_B; ++i) {
                align_B[target_segment1 * COL_B + i] += source[i];
                align_B[target_segment2 * COL_B + i] += source[i];
            }
        }

        // for (int segment = tile_offset; segment < tile_offset + wl[tidx].num_Merge; ++segment) {
        //     if (merging[segment].dst & 0x80000000) {
        //         DATA_TYPE target = merging[segment].dst & 0x0000000f;
        //         memset(merge_split + (target * COL_B), 0, COL_B * sizeof(float));
        //     }
        // }
        end_Time(TIME_TEST2);
    }

    end_Time(TIME_ALIGN);

    return padded_segment;
}


struct thread_args {
    VAL_TYPE *align_A;
    VAL_TYPE *align_B;
    VAL_TYPE *C;
    int padded_segment;
    int binded_tile;
    int mode;
    // int tile_index;
    // uint8_t status;
    // uint8_t split;
    // uint16_t merge;
};
typedef struct thread_args thread_args;

void *manage_pim(void *args) {
    thread_args *worker_args = (thread_args *)args;
    VAL_TYPE *align_A = worker_args->align_A;
    VAL_TYPE *align_B = worker_args->align_B;
    VAL_TYPE *C = worker_args->C;
    int mode = worker_args->mode;
    uint32_t A_ROW_PAD = tmat.nr * REG_SIZE; // M, p
    uint32_t A_COL_PAD = (mat.nc + REG_SIZE - 1) / REG_SIZE * REG_SIZE; // K, q

    start_Time(TIME_OVERLAP);
    if (check_option(mode, BINDING)) {
        int cur_tile = 0;
        int alter = 0;
        uint32_t *num_segments = (uint32_t*)malloc(tmat.nr << 2);
        while (cur_tile < tmat.nr) {
            VAL_TYPE *cur_align_A = align_A + (NUM_BANK * A_COL_PAD * alter);
            VAL_TYPE *cur_align_B = align_B + (A_COL_PAD * COL_B * alter);

            while (status != WORKING) usleep(1);
            int padded_segment = worker_args->padded_segment;
            uint32_t num_tile = worker_args->binded_tile;
            for (int tile = 0; tile < num_tile; ++tile) {
                if (check_option(mode, MERGING)) num_segments[tile] = wl[cur_tile + tile].num_Segment;
                else num_segments[tile] = tmat.ptr[cur_tile + tile + 1] - tmat.ptr[cur_tile + tile];
            }

            status = WAITING;
            // binded_matmul_cpu(cur_align_A, cur_align_B, C + (cur_tile * COL_B * NUM_BANK), num_segments, num_tile);
            emulate_gemm_pim(align_A, align_B, C + (cur_tile * COL_B * NUM_BANK), padded_segment);
            cur_tile += num_tile;
            alter = 1 - alter;
        }
        free(num_segments);
    }
    else {
        for (int tile = 0; tile < tmat.nr; ++tile) {
            // int num_segment = tmat.ptr[tile + 1] - tmat.ptr[tile];
            // if (check_option(mode, MERGING)) num_segment = wl[tile].num_Segment;
            // int num_gemm = (num_segment + NUM_BANK - 1) / NUM_BANK; 
            // int padded_segment = (num_segment + NUM_BANK - 1) / NUM_BANK * NUM_BANK;

            VAL_TYPE *cur_align_A = align_A + (NUM_BANK * A_COL_PAD * (tile & 0x1));
            VAL_TYPE *cur_align_B = align_B + (A_COL_PAD * COL_B * (tile & 0x1));

            while (status != WORKING) usleep(1);
            int padded_segment = worker_args->padded_segment;

            status = WAITING;
            // matmul_cpu(cur_align_A, cur_align_B, C + (tile * COL_B * NUM_BANK), padded_segment);
            emulate_gemm_pim(align_A, align_B, C + (tile * COL_B * NUM_BANK), padded_segment);
        }
    }
    end_Time(TIME_OVERLAP);

    return NULL;
}

void tiled_spmm_with_tcsr(VAL_TYPE *align_A, VAL_TYPE *B, VAL_TYPE *align_B, VAL_TYPE *C, uint32_t M, uint32_t K, uint32_t N, int mode) {
    init_SpMM(mode);

    if (check_option(mode, OVERLAPPING)) {
        // uint32_t A_ROW_PAD = (mat.nr + NUM_BANK - 1) / NUM_BANK * NUM_BANK;
        uint32_t A_COL_PAD = (mat.nc + NUM_BANK - 1) / NUM_BANK * NUM_BANK; // K, q
        thread_args args = {align_A, align_B, C, 0, 0, mode};
        status = WAITING;
        pthread_t pim_worker;
        pthread_create(&pim_worker, NULL, manage_pim, &args);
        if (check_option(mode, BINDING)) {
            int tile = 0, binded_tile = 0, alter = 0;
            while (tile < tmat.nr) {
                // printf("Current Tile: %d / %d\n", tile, tmat.nr);
                while (status != WAITING) usleep(1);
        
                int num_tile = 0;
                int acc_segment = 0;
                while (acc_segment < BIND_CRITERIA & (tile + num_tile < tmat.nr)) {
                    if (check_option(mode, MERGING)) {
                        merge_Tile(wl, merging, splitting, tile + num_tile);
                    }
                    if (check_option(mode, PARALLEL)) 
                        acc_segment += align_tile_parallel(align_A + (acc_segment * NUM_BANK), B, align_B + (acc_segment * COL_B), tile + num_tile, mode);
                    else
                        acc_segment += align_tile(align_A + (acc_segment * NUM_BANK), B, align_B + (acc_segment * COL_B), tile + num_tile, mode);
                    num_tile++;
                }
                
                // args.padded_segment = align_tile_binded(align_A + (NUM_BANK * A_COL_PAD * alter), pim_A + (NUM_BANK * A_COL_PAD * alter), B, align_B + (A_COL_PAD * COL_B * alter), pim_B + (A_COL_PAD * COL_B * alter), tile, tile + num_tile, mode);
                args.binded_tile = num_tile;
                args.padded_segment = acc_segment;
                status = WORKING;
                tile += num_tile;
                binded_tile++;
                alter = 1 - alter;
            }
            print_Binding(binded_tile);
        }
        else {
            for (int tile = 0; tile < tmat.nr; ++tile) {
                // printf("Current Tile: %d / %d\n", tile, tmat.nr);
                while (status != WAITING) usleep(1);
        
                if (check_option(mode, MERGING)) merge_Tile(wl, merging, splitting, tile);
                args.padded_segment = align_tile(align_A + (NUM_BANK * A_COL_PAD * (tile & 0x1)), B, align_B + (A_COL_PAD * COL_B * (tile & 0x1)), tile, mode);
                status = WORKING;
            }
        }
        pthread_join(pim_worker, NULL);
    }
    else {
        if (check_option(mode, BINDING)) {
            int tile = 0, binded_tile = 0;
            uint32_t *num_segments = (uint32_t*)malloc(tmat.nr << 2);
            while (tile < tmat.nr) {
        
                int num_tile = 0;
                int acc_segment = 0;
                while (acc_segment < BIND_CRITERIA & (tile + num_tile < tmat.nr)) {
                    if (check_option(mode, MERGING)) {
                        merge_Tile(wl, merging, splitting, tile + num_tile);
                    }
                    num_segments[num_tile] = align_tile(align_A + (acc_segment * NUM_BANK), B, align_B + (acc_segment * COL_B), tile + num_tile, mode);
                    acc_segment += num_segments[num_tile];
                    num_tile++;
                }
                
                // uint32_t padded_segment = align_tile_binded(align_A, pim_A, B, align_B, pim_B, tile, tile + num_tile, mode);
                // binded_matmul_cpu(align_A, align_B, C + (tile * N * NUM_BANK), num_segments, num_tile);
                emulate_gemm_pim(align_A, align_B, C + (tile * N * NUM_BANK), acc_segment);
                tile += num_tile;
                binded_tile++;
            }
            free(num_segments);
            print_Binding(binded_tile);
        }
        else {
            for (int tile = 0; tile < tmat.nr; ++tile) {
                // printf("Current Tile: %d / %d\n", tile, tmat.nr);
                if (check_option(mode, MERGING)) merge_Tile(wl, merging, splitting, tile);
                int padded_segment = align_tile(align_A, B, align_B, tile, mode);
                // int padded_segment = align_tile_parallel(align_A, B, align_B, tile, mode);

                // printf("TILE %d: %d\n", tile, padded_segment);
                matmul_cpu(align_A, align_B, C, padded_segment);
                // matmul_cpu(align_A, align_B, C + (tile * N * NUM_BANK), padded_segment);
                // emulate_gemm_pim(align_A, align_B, C + (tile * N * NUM_BANK), padded_segment);
            }
        }
    }
    
    end_SpMM(mode);
}



void tiled_spmm(VAL_TYPE *align_A, VAL_TYPE *B, VAL_TYPE *align_B, VAL_TYPE *C, uint32_t M, uint32_t K, uint32_t N, int mode) {
    // Tiled SpMM with tCSR Format
    if (check_option(mode, TCSR_SPMM)) tiled_spmm_with_tcsr(align_A, B, align_B, C, M, K, N, mode);
}
