#include "SpMM.hpp"

void move_Rows();
void reorder_Rows();
void reorder_Rows2();
void check_Segments();
void check_Segments_parallel();
void count_Segments();

void make_TiledCSR(int mode) {
    move_Rows();
    if (check_option(mode, TCSR_SPMM)) {
        // if (check_option(mode, PARALLEL)) check_Segments_parallel();
        // else check_Segments();
        check_Segments_parallel();
        // check_Segments();

    }
    else if (check_option(mode, CSR_SPMM)) count_Segments();
    // check_Segments();
}

void move_Rows() {
    int first_offset = 0;
    int last_offset = mat.nr - 1;
    DATA_TYPE* temp_nnz = (DATA_TYPE*)malloc(mat.nr << SHIFT_BIT);


    for (int i = 0; i < mat.nr; ++i) {
        if (mat.nnz[i] == 0) {
            temp_nnz[last_offset] = 0;
            last_offset--;
        }
        else {
            temp_nnz[first_offset] = mat.nnz[i];
            first_offset++;
        }
    }
    for (int i = 0; i < mat.nr; ++i) {
        mat.ptr[i + 1] = mat.ptr[i] + temp_nnz[i];
    }
    free(mat.nnz);
    mat.nnz = temp_nnz;

    tmat.nr = (first_offset + NUM_BANK - 1) >> 4;


    int original_nr = mat.nr;
    mat.nr = first_offset;
    puts(" --- Remove Empty Rows --------------------------------------------");
    printf("| Number of Row: %d -> %d (%.2f%%)\n\n", original_nr, first_offset, (float)100 * first_offset / original_nr);
}

void count_Segments() {
    tmat.ptr = (DATA_TYPE*)malloc((tmat.nr + 1) * sizeof(DATA_TYPE));
    // tmat.mapping = (DATA_TYPE*)malloc(tmat.nr * mat.nc << SHIFT_BIT);
    // tmat.num = (uint8_t*)calloc(tmat.nr * mat.nc, 1);
    uint8_t *exist = (uint8_t*)malloc(mat.nc);
    memset(exist, 0x00, mat.nc);
    int num_GEMM = 0;

    start_Time(TIME_TCSR);
    tmat.ptr[0] = 0;
    for (int tile = 0; tile < tmat.nr; ++tile) {
        int target = tile < tmat.nr - 1 ? (tile << 4) + 16 : mat.nr;
        int count = 0;
        for (int element = mat.ptr[tile << 4]; element < mat.ptr[target]; ++element) {
            if (!exist[mat.idx[element]]) {
                count++;
                exist[mat.idx[element]] = 1;
            }
        }
        // tmat.num[tile] = count;
        tmat.ptr[tile + 1] = tmat.ptr[tile] + count;
        num_GEMM += (count + NUM_BANK - 1) >> 4;
        
        for (int element = mat.ptr[tile << 4]; element < mat.ptr[target]; ++element) {
            exist[mat.idx[element]] = 0;
        }
    }
    tmat.ne = tmat.ptr[tmat.nr];
    
    end_Time(TIME_TCSR);
    free(exist);

    puts(" --- Count Number of Segments -------------------------------------");
    printf("| Number of Segments: %d\n", tmat.ne);
    printf("| Number of GEMM Units: %d\n\n", num_GEMM);
}

void check_Segments_parallel() {
#define NUM_THREAD 16
    tmat.ptr = (DATA_TYPE*)malloc(tmat.nr + 1 << SHIFT_BIT);
    tmat.mapping = (DATA_TYPE*)malloc(tmat.nr * mat.nc << SHIFT_BIT);
    tmat.num = (uint8_t*)calloc(tmat.nr * mat.nc, 1);
    tmat.row_hash = (uint16_t*)calloc(tmat.nr * mat.nc, 2);

    uint8_t *exist = (uint8_t*)malloc(mat.nc * NUM_THREAD);
    memset(exist, 0x00, mat.nc * NUM_THREAD);
    DATA_TYPE *nnzs = (DATA_TYPE*)calloc(tmat.nr << SHIFT_BIT, 1);

    DATA_TYPE *mapping = (DATA_TYPE*)malloc(mat.nc * NUM_THREAD << SHIFT_BIT);
    memset(mapping, 0xff, mat.nc * NUM_THREAD << SHIFT_BIT);

    int num_GEMM = 0;
    start_Time(TIME_TCSR);
#pragma omp parallel for num_threads(NUM_THREAD) schedule(guided)
    for (int tile = 0; tile < tmat.nr; ++tile) {
        uint8_t *my_exist = exist + (omp_get_thread_num() * mat.nc);
        int target = tile < tmat.nr - 1 ? (tile << 4) + 16 : mat.nr;
        int count = 0;
        for (int element = mat.ptr[tile << 4]; element < mat.ptr[target]; ++element) {
            if (!my_exist[mat.idx[element]]) {
                count++;
                my_exist[mat.idx[element]] = 1;
            }
        }
        nnzs[tile] = count;
        
        for (int element = mat.ptr[tile << 4]; element < mat.ptr[target]; ++element) {
            my_exist[mat.idx[element]] = 0;
        }
    }
    tmat.ptr[0] = 0;
    for (int tile = 0; tile < tmat.nr; ++tile) {
        num_GEMM += (nnzs[tile] + NUM_THREAD - 1) >> 4;
        tmat.ptr[tile + 1] = tmat.ptr[tile] + nnzs[tile];
    } 
    tmat.ne = tmat.ptr[tmat.nr];

#pragma omp parallel for num_threads(NUM_THREAD) schedule(guided)
    for (int tile = 0; tile < tmat.nr; ++tile) {
        DATA_TYPE *my_mapping = mapping + (omp_get_thread_num() * mat.nc);
        DATA_TYPE first = 0;
        int target = tile < tmat.nr - 1 ? (tile << 4) + 16 : mat.nr;
        for (int row = tile << 4; row < target; ++row) {
            uint16_t hash_val = 0x1 << ((uint16_t)row & 0x000f);
            for (int element = mat.ptr[row]; element < mat.ptr[row + 1]; ++element) {
                DATA_TYPE cur_col = mat.idx[element];
                if (my_mapping[cur_col] == 0xffffffff) {
                    my_mapping[cur_col] = first;
                    tmat.mapping[tmat.ptr[tile] + first] = cur_col;
                    first++;
                }
                int mapped_seg = tmat.ptr[tile] + my_mapping[cur_col];
                tmat.row_hash[mapped_seg] += hash_val;
                tmat.num[mapped_seg]++;
            }
        }
        tmat.ptr[tile + 1] = tmat.ptr[tile] + first;

        for (int i = 0; i < first; ++i) my_mapping[tmat.mapping[tmat.ptr[tile] + i]] = 0xffffffff;
    }

    tmat.ne = tmat.ptr[tmat.nr];

    end_Time(TIME_TCSR);
    free(mapping);
    free(exist);
    free(nnzs);

// Check integrity
    // DATA_TYPE count = 0;
    // for (int i = 0; i < tmat.ne; ++i) {
    //     count += tmat.num[i];
    // }
    // assert(mat.ne == count);
	// printf("%d vs %d\n", mat.ne, count);

    puts(" --- Build Tiled CSR Format ---------------------------------------");
    printf("| Number of Segments: %d\n", tmat.ne);
    printf("| Number of GEMM Units: %d\n\n", num_GEMM);
}

void check_Segments() {
    tmat.ptr = (DATA_TYPE*)malloc(tmat.nr + 1 << SHIFT_BIT);
    tmat.mapping = (DATA_TYPE*)malloc(tmat.nr * mat.nc << SHIFT_BIT);
    tmat.num = (uint8_t*)calloc(tmat.nr * mat.nc, 1);
    tmat.row_hash = (uint16_t*)calloc(tmat.nr * mat.nc, 2);
    DATA_TYPE* mapping = (DATA_TYPE*)malloc(mat.nc << SHIFT_BIT);
    memset(mapping, 0xff, mat.nc << SHIFT_BIT);

    int num_GEMM = 0;
    start_Time(TIME_TCSR);

    tmat.ptr[0] = 0;
    for (int tile = 0; tile < tmat.nr; ++tile) {
        DATA_TYPE first = 0;
        int target = tile < tmat.nr - 1 ? (tile << 4) + 16 : mat.nr;
        for (int row = tile << 4; row < target; ++row) {
            uint16_t hash_val = 0x1 << ((uint16_t)row & 0x000f);
            for (int element = mat.ptr[row]; element < mat.ptr[row + 1]; ++element) {
                DATA_TYPE cur_col = mat.idx[element];
                if (mapping[cur_col] == 0xffffffff) {
                    mapping[cur_col] = first;
                    tmat.mapping[tmat.ptr[tile] + first] = cur_col;
                    first++;
                }
                int mapped_seg = tmat.ptr[tile] + mapping[cur_col];
                tmat.row_hash[mapped_seg] += hash_val;
                tmat.num[mapped_seg]++;
            }
        }
        tmat.ptr[tile + 1] = tmat.ptr[tile] + first;
        num_GEMM += (first + NUM_BANK - 1) >> 4;

        for (int i = 0; i < first; ++i) mapping[tmat.mapping[tmat.ptr[tile] + i]] = 0xffffffff;
    }

    tmat.ne = tmat.ptr[tmat.nr];

    end_Time(TIME_TCSR);
    free(mapping);

// Check integrity
    DATA_TYPE count = 0;
    for (int i = 0; i < tmat.ne; ++i) {
        count += tmat.num[i];
    }
    // assert(mat.ne == count);
	printf("%d vs %d\n", mat.ne, count);

    puts(" --- Build Tiled CSR Format ---------------------------------------");
    printf("| Number of Segments: %d\n", tmat.ne);
    printf("| Number of GEMM Units: %d\n\n", num_GEMM);
}

inline void swap(DATA_TYPE *array, int a, int b, int num) {
    for (int i = 0; i < num; ++i) {
        DATA_TYPE temp = array[a + i];
        array[a + i] = temp;
        array[a + i] = array[b + i];
        array[b + i] = temp;
    }
    // printf("SWAP: %d (%d'th row in Tile %d) and %d with %d size\n", a, a & 0xf, a >> 4, b, num);
}

#define CRITERIA (float)3.0f

// void reorder_Rows2() {
//     char *dense = (char*)calloc(mat.nc, 1);
//     DATA_TYPE *group_ptr = (DATA_TYPE*)malloc(mat.nr + 1 << SHIFT_BIT);
//     DATA_TYPE *group_idx = (DATA_TYPE*)malloc(mat.nr << SHIFT_BIT);
//     mat.mapping = (DATA_TYPE*)malloc(mat.nr << SHIFT_BIT);
//     for (int i = 0; i < mat.nr; ++i) mat.mapping[i] = i;

//     start_Time(TIME_REORDER);
//     start_Time(TIME_TEST);
//     int row = 0, num_group = 0;
//     group_ptr[0] = 0;
//     while (row < mat.nr) {
//         DATA_TYPE total_score = mat.nnz[row];
//         group_idx[row] = num_group;
//         group_ptr[num_group] = row;
//         for (int element = mat.ptr[row]; element < mat.ptr[row + 1]; ++element) {
//             dense[mat.idx[element]]++;
//         }
//         row++;

//         // bool search = true;
//         bool search = row < mat.nr ? true : false;
//         int level = 1;
//         while (search) {
//             DATA_TYPE cur_score = 0;
//             for (int element = mat.ptr[row]; element < mat.ptr[row + 1]; ++element) {
//                 cur_score += dense[mat.idx[element]];
//             }
//             // if (cur_score * CRITERIA >= total_score) {
//             if (floor(cur_score * CRITERIA) >= total_score) {
//                 for (int element = mat.ptr[row]; element < mat.ptr[row + 1]; ++element) {
//                     dense[mat.idx[element]]++;
//                 }
//                 total_score += mat.nnz[row];
//                 group_idx[row] = num_group;
//                 level++;
//                 row++;
//                 if (row == mat.nr || level == NUM_BANK) {
//                     for (int element = mat.ptr[row - level]; element < mat.ptr[row]; ++element) {
//                         dense[mat.idx[element]] = 0;
//                     }
//                     search = false;
//                 }
//             }
//             else {
//                 num_group++;
//                 search = false;
//                 for (int element = mat.ptr[row - level]; element < mat.ptr[row]; ++element) {
//                     dense[mat.idx[element]] = 0;
//                 }
//             }
//         }
//     }
//     end_Time(TIME_TEST);

//     start_Time(TIME_TEST2);
//     int total_swap = 0;
//     for (int group = 0; group < num_group - 1; ++group) {
//         int first_row_in_cur_group = group_ptr[group];
//         int last_row_in_cur_group = group_ptr[group + 1] - 1;

//         if ((first_row_in_cur_group & 0xf) > (last_row_in_cur_group & 0xf)) {
//             int num_row_in_next_group = (last_row_in_cur_group & 0xf) + 1;
//             int first_group_in_cur_tile = group_idx[first_row_in_cur_group & 0xfffffff0];
//             first_group_in_cur_tile = group_ptr[first_group_in_cur_tile] < first_row_in_cur_group ? first_group_in_cur_tile + 1 : first_group_in_cur_tile;
//             int left = first_group_in_cur_tile;

//             for (int i = first_group_in_cur_tile; i < group; ++i) {
//                 int group_size = group_ptr[i + 1] - group_ptr[left];
//                 while (group_size > num_row_in_next_group) {
//                     left++;
//                     group_size = group_ptr[i + 1] - group_ptr[left];
//                 }
//                 if (group_size == num_row_in_next_group) {
//                     swap(mat.mapping, group_ptr[left], last_row_in_cur_group & 0xfffffff0, num_row_in_next_group);
//                     total_swap += num_row_in_next_group;
//                     break;
//                 }
//                 // else {
//                 //     while (group_size > num_row_in_next_group) {
//                 //         left++;
//                 //         group_size = group_ptr[i + 1] - group_ptr[left];
//                 //     }
//                 //     if (group_size == num_row_in_next_group) {
//                 //         swap(mat.mapping, group_ptr[left], last_row_in_cur_group & 0xfffffff0, num_row_in_next_group);
//                 //         total_swap += num_row_in_next_group;
//                 //         break;
//                 //     }
//                 // }
//             }

//         }

//         // int cur_group = group_ptr[group];
//         // int next_group = group_ptr[group + 1];

//         // if ((cur_group & 0x0000000f) > (next_group & 0x0000000f)) {
//         //     // int cur_num = NUM_BANK - (cur_group & 0x0000000f);
//         //     int next_num = (next_group & 0x0000000f) + 1;
//         //     printf("%d ~ %d\n", cur_group, next_group - 1);
            
//         //     // Pick rows in upper tile
//         //     int start_group = group_idx[cur_group & 0xfffffff0];
//         //     int left = group_ptr[start_group] < cur_group ? start_group + 1 : start_group;
//         //     // int left = start_group;

//         //     for (int i = start_group; i < cur_group; ++i) {
//         //         int group_size = group_ptr[i + 1] - group_ptr[left];
//         //         if (group_size == next_num) {
//         //             swap(mat.mapping, group_ptr[left], next_group & 0xfffffff0, next_num);
//         //             total_swap += next_num;
//         //             break;
//         //         }
//         //         else {
//         //             while (group_size > next_num) {
//         //                 left++;
//         //                 group_size = group_ptr[i + 1] - group_ptr[left];
//         //             }
//         //             if (group_size == next_num) {
//         //                 swap(mat.mapping, group_ptr[left], next_group & 0xfffffff0, next_num);
//         //                 total_swap += next_num;
//         //                 break;
//         //             }
//         //         }
//         //     }
//         // }
//     }

//     end_Time(TIME_TEST2);
//     end_Time(TIME_REORDER);

//     free(dense);
//     free(group_ptr);
//     free(group_idx);

//     puts(" --- Reorder Rows with Grouping -----------------------------------");
//     printf("| Number of Groups: %d\n", num_group);
//     printf("| Number of Exchanges: %d\n\n", total_swap);
// }

// struct tile_info {
//     int target_idx;
//     int target_score;
//     int8_t *dense;
// };

// void reorder_Rows() {
//     mat.mapping = (DATA_TYPE*)malloc(mat.nr << SHIFT_BIT);
//     for (int i = 0; i < mat.nr; ++i) mat.mapping[i] = i;

//     struct tile_info tinfo[2];
//     for (int i = 0; i < 2; ++i) tinfo[i].dense = (int8_t*)calloc(mat.nc, 1);

//     start_Time(TIME_REORDER);

//     // Initialize dense array with range of tiles
//     for (int element = mat.ptr[0]; element < mat.ptr[NUM_BANK]; ++element) {
//         tinfo[0].dense[mat.idx[element]]++;
//     }

//     int total_count = 0;

//     // Exchange rows which exists in two adjacent tiles
//     tile_info *cur_tile = tinfo + 1, *next_tile = tinfo;
//     for (int tile = 0; tile < tmat.nr - 1; ++tile) {
//         tile_info *temp = cur_tile;
//         cur_tile = next_tile;
//         next_tile = temp;
//         int range = ((tile + 2) << 4) < mat.nr ? (tile + 2) << 4 : mat.nr;

//         // Make dense array of next tile
//         memset(next_tile->dense, 0, mat.nc);
//         for (int element = mat.ptr[(tile + 1) << 4]; element < mat.ptr[range]; ++element) {
//             next_tile->dense[mat.idx[element]]++;
//         }

//         // Find a row which creates the most segments by itself (in current tile)
//         cur_tile->target_score = 0x80000001;
//         cur_tile->target_idx = 0xffffffff;
//         for (int row = tile << 4; row < (tile + 1) << 4; ++row) {
//             int score = 0;    // The number of elements which create segments by themselves
//             int mapped_row = mat.mapping[row];
//             for (int element = mat.ptr[mapped_row]; element < mat.ptr[mapped_row + 1]; ++element) {
//                 if (cur_tile->dense[mat.idx[element]] == 1) score++;
//                 if (next_tile->dense[mat.idx[element]] == 0) score--;
//             }
//             if (score > cur_tile->target_score) {
//                 cur_tile->target_score = score;
//                 cur_tile->target_idx = row;
//             }
//         }
//         for (int element = mat.ptr[mat.mapping[cur_tile->target_idx]]; element < mat.ptr[mat.mapping[cur_tile->target_idx] + 1]; ++element) {
//             cur_tile->dense[mat.idx[element]]--;
//             // if (cur_tile->dense[mat.idx[element] == 1]) cur_tile->dense[mat.idx[element]] = 0x40;
//         }

//         // Find a row which creates the most segments by itself (in next tile)
//         next_tile->target_score = 0x80000001;
//         next_tile->target_idx = 0xffffffff;
//         int selected_row = mat.mapping[cur_tile->target_idx];
//         for (int row = (tile + 1) << 4; row < range; ++row) {
//             int score = 0;    // The number of elements which create segments by themselves
//             for (int element = mat.ptr[row]; element < mat.ptr[row + 1]; ++element) {
//                 if (next_tile->dense[mat.idx[element]] == 1) score++;
//                 if (cur_tile->dense[mat.idx[element]] == 0) score--;
//                 // if (next_tile->dense[mat.idx[element]] == 1) score++;
//                 // if (cur_tile->dense[mat.idx[element]] == 0 || cur_tile->dense[mat.idx[element]] == 0x40) score--;
//                 // if (cur_tile->dense[mat.idx[element]] == 0x40 && next_tile->dense[mat.idx[element] == 1]) score--;
                
//             }
//             // score -= compare_row(mat.mapping[cur_tile->target_idx], row);
//             if (score > next_tile->target_score) {
//                 next_tile->target_score = score;
//                 next_tile->target_idx = row;
//             }
//         }
//         for (int element = mat.ptr[mat.mapping[next_tile->target_idx]]; element < mat.ptr[mat.mapping[next_tile->target_idx] + 1]; ++element) {
//             next_tile->dense[mat.idx[element]]--;
//         }
//         cur_tile->target_score = 0;
//         for (int element = mat.ptr[mat.mapping[cur_tile->target_idx]]; element < mat.ptr[mat.mapping[cur_tile->target_idx] + 1]; ++element) {
//             if (cur_tile->dense[mat.idx[element]] == 1) cur_tile->target_score++;
//             if (next_tile->dense[mat.idx[element]] == 0) cur_tile->target_score--;
//         }

//         // Exchange if the number of segments decreases
//         if (cur_tile->target_score > 0 && next_tile->target_score > 0) {
//             swap(mat.mapping, cur_tile->target_idx, next_tile->target_idx, 1);
//             DATA_TYPE cur_row = mat.mapping[next_tile->target_idx];
//             DATA_TYPE next_row = mat.mapping[cur_tile->target_idx];
//             for (int element = mat.ptr[cur_row]; element < mat.ptr[cur_row + 1]; ++element) next_tile->dense[mat.idx[element]]++;
//             // for (int element = mat.ptr[next_row]; element < mat.ptr[next_row + 1]; ++element) next_tile->dense[mat.idx[element]]--;

//             // for (int j = 0; j < mat.nnz[cur_row]; ++j) printf("\t(%d %d)\n", mat.idx[mat.ptr[cur_row] + j], mat.idx[mat.ptr[next_row] + j]);
//             // total_score += cur_tile->target_score + next_tile->target_score;
//             total_count++;
//         }
//         else {
//             for (int element = mat.ptr[mat.mapping[next_tile->target_idx]]; element < mat.ptr[mat.mapping[next_tile->target_idx] + 1]; ++element)
//                 next_tile->dense[mat.idx[element]]++;
//         }

//     }
//     // printf("Total score: %d\n", total_score);

//     end_Time(TIME_REORDER);
//     for (int i = 0; i < 2; ++i) free(tinfo[i].dense);

//     puts(" --- Reorder Rows Between Adjacent Tiles---------------------------");
//     printf("| Number of Exchanges: %d\n\n", total_count);
// }


// void reorder_Rows() {
//     mat.mapping = (DATA_TYPE*)malloc(mat.nr << SHIFT_BIT);
//     for (int i = 0; i < mat.nr; ++i) mat.mapping[i] = i;

//     start_Time(TIME_REORDER);

//     struct tile_info tinfo[2];
//     for (int i = 0; i < 2; ++i) tinfo[i].dense = (int8_t*)calloc(mat.nc, 1);

//     // Initialize dense array with range of tiles
//     for (int element = mat.ptr[0]; element < mat.ptr[NUM_BANK]; ++element) {
//         tinfo[0].dense[mat.idx[element]]++;
//     }

//     int total_score = 0;

//     // Exchange rows which exists in two adjacent tiles
//     tile_info *cur_tile = tinfo + 1, *next_tile = tinfo;
//     for (int tile = 1; tile < tmat.nr - 1; ++tile) {
//         tile_info *temp = cur_tile;
//         cur_tile = next_tile;
//         next_tile = temp;
//         int range = ((tile + 2) << 4) < mat.nr ? (tile + 2) << 4 : mat.nr;

//         // Make dense array of next tile
//         memset(next_tile->dense, 0, mat.nc);
//         for (int element = mat.ptr[(tile + 1) << 4]; element < mat.ptr[range]; ++element) {
//             next_tile->dense[mat.idx[element]]++;
//         }

//         int changed = 1;
//         // while (changed) {
//         for (int i = 0; i < NUM_BANK & changed; ++i) {
//             changed = 0;

//             // Find a row which creates the most segments by itself (in current tile)
//             cur_tile->target_score = 0x80000001;
//             cur_tile->target_idx = 0xffffffff;
//             for (int row = tile << 4; row < (tile + 1) << 4; ++row) {
//                 int score = 0;    // The number of elements which create segments by themselves
//                 int mapped_row = mat.mapping[row];
//                 for (int element = mat.ptr[mapped_row]; element < mat.ptr[mapped_row + 1]; ++element) {
//                     if (cur_tile->dense[mat.idx[element]] == 1) score++;
//                     if (next_tile->dense[mat.idx[element]] == 0) score--;
//                 }
//                 if (score > cur_tile->target_score) {
//                     cur_tile->target_score = score;
//                     cur_tile->target_idx = row;
//                 }
//             }
//             for (int element = mat.ptr[mat.mapping[cur_tile->target_idx]]; element < mat.ptr[mat.mapping[cur_tile->target_idx] + 1]; ++element) {
//                 cur_tile->dense[mat.idx[element]]--;
//             }

//             // Find a row which creates the most segments by itself (in next tile)
//             next_tile->target_score = 0x80000001;
//             next_tile->target_idx = 0xffffffff;
//             for (int row = (tile + 1) << 4; row < range; ++row) {
//                 int score = 0;    // The number of elements which create segments by themselves
//                 int mapped_row = mat.mapping[row];
//                 for (int element = mat.ptr[mapped_row]; element < mat.ptr[mapped_row + 1]; ++element) {
//                     if (cur_tile->dense[mat.idx[element]] == 0) score--;
//                     if (next_tile->dense[mat.idx[element]] == 1) score++;
//                 }
//                 if (score > next_tile->target_score) {
//                     next_tile->target_score = score;
//                     next_tile->target_idx = row;
//                 }
//             }

//             // Exchange if the number of segments decreases
//             if (cur_tile->target_score > 0 && next_tile->target_score > 0) {
//             // if (cur_tile->target_score + next_tile->target_score > 0) {
//                 swap(mat.mapping, cur_tile->target_idx, next_tile->target_idx);
//                 DATA_TYPE cur_row = mat.mapping[cur_tile->target_idx];
//                 DATA_TYPE next_row = mat.mapping[next_tile->target_idx];
//                 for (int element = mat.ptr[cur_row]; element < mat.ptr[cur_row + 1]; ++element) {
//                     // cur_tile->dense[mat.idx[element]]--;
//                     next_tile->dense[mat.idx[element]]++;
//                 }

//                 for (int element = mat.ptr[next_row]; element < mat.ptr[next_row + 1]; ++element) {
//                     cur_tile->dense[mat.idx[element]]++;
//                     next_tile->dense[mat.idx[element]]--;
//                 }
//                 printf("Exchanged %d and %d (score: %d and %d)\n", cur_tile->target_idx, next_tile->target_idx, cur_tile->target_score, next_tile->target_score);
//                 total_score += cur_tile->target_score + next_tile->target_score;
//                 changed = 1;
//             }
//         }

//     }
//     printf("Total score: %d\n", total_score);

//     end_Time(TIME_REORDER);
//     for (int i = 0; i < 2; ++i) free(tinfo[i].dense);
// }
