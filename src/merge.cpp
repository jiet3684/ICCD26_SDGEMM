#include "SpMM.hpp"

uint16_t *hash_values;
uint8_t *merged;
DATA_TYPE *col_src;

DATA_TYPE *first;
DATA_TYPE *mapping;
char *splitted;
#define SEGMENT_SPLIT true;

void init_Merge(int mode) {
    first = (DATA_TYPE*)malloc(0xffff << SHIFT_BIT);
    memset(first, 0xff, 0xffff << SHIFT_BIT);
    mapping = (DATA_TYPE*)malloc(tmat.ne << SHIFT_BIT);
    splitted = (char*)malloc(tmat.ne);
    memset(splitted, 0xff, tmat.ne);
}

void end_Merge(int mode) {
    free(first);
    free(splitted);
    free(mapping);
}


inline bool split_Segments(DATA_TYPE *first, split *splitting, DATA_TYPE segment, DATA_TYPE num) {
    if (num > 4 || num <= 1) return false;

    uint16_t hash_val = tmat.row_hash[segment];
    if (num == 2) {
        uint32_t idx1, idx2;
        int i;
        for (i = 0; i < NUM_BANK; ++i) {
            uint16_t temp_val = hash_val & (0x0001 << i);
            if (temp_val != 0x0000) {
                idx1 = temp_val;
                break;
            }
        }
        if (i >= NUM_BANK) return false;
        for (i = idx1 + 1; i < NUM_BANK; ++i) {
            uint16_t temp_val = hash_val & (0x0001 << i);
            if (temp_val != 0x0000) {
                idx2 = temp_val;
                break;
            }
        }
        if (i >= NUM_BANK) return false;

        if (first[idx1] != 0xffffffff && first[idx2] != 0xffffffff) {
            splitting->dst1 = first[idx1];
            splitting->dst2 = first[idx2];
        }
        else return false;

    }
else if (num == 3) {
        uint32_t idx1, idx2, idx3;
        int i;
        for (i = 0; i < NUM_BANK; ++i) {
            uint16_t temp_val = hash_val & (0x0001 << i);
            if (temp_val != 0x0000) {
                idx1 = temp_val;
                break;
            }
        }
        if (i >= NUM_BANK) return false;
        for (i = idx1 + 1; i < NUM_BANK; ++i) {
            uint16_t temp_val = hash_val & (0x0001 << i);
            if (temp_val != 0x0000) {
                idx2 = temp_val;
                break;
            }
        }
        if (i >= NUM_BANK) return false;
        for (i = idx2 + 1; i < NUM_BANK; ++i) {
            uint16_t temp_val = hash_val & (0x0001 << i);
            if (temp_val != 0x0000) {
                idx3 = temp_val;
                break;
            }
        }
        if (i >= NUM_BANK) return false;

        if (first[idx1 + idx2] != 0xffffffff && first[idx3] != 0xffffffff) {
            splitting->dst1 = first[idx1 + idx2];
            splitting->dst2 = first[idx3];
        }
        else if (first[idx1 + idx3] != 0xffffffff && first[idx2] != 0xffffffff) {
            splitting->dst1 = first[idx1 + idx3];
            splitting->dst2 = first[idx2];
        }
        else if (first[idx2 + idx3] != 0xffffffff && first[idx1] != 0xffffffff) {
            splitting->dst1 = first[idx2 + idx3];
            splitting->dst2 = first[idx1];
        }
        else return false;

    }
    else if (num == 4) {
        uint32_t idx1, idx2, idx3, idx4;
        int i;
        for (i = 0; i < NUM_BANK; ++i) {
            uint16_t temp_val = hash_val & (0x0001 << i);
            if (temp_val != 0x0000) {
                idx1 = temp_val;
                break;
            }
        }
        if (i >= NUM_BANK) return false;
        for (i = idx1 + 1; i < NUM_BANK; ++i) {
            uint16_t temp_val = hash_val & (0x0001 << i);
            if (temp_val != 0x0000) {
                idx2 = temp_val;
                break;
            }
        }
        if (i >= NUM_BANK) return false;
        for (i = idx2 + 1; i < NUM_BANK; ++i) {
            uint16_t temp_val = hash_val & (0x0001 << i);
            if (temp_val != 0x0000) {
                idx3 = temp_val;
                break;
            }
        }
        if (i >= NUM_BANK) return false;
        for (i = idx3 + 1; i < NUM_BANK; ++i) {
            uint16_t temp_val = hash_val & (0x0001 << i);
            if (temp_val != 0x0000) {
                idx4 = temp_val;
                break;
            }
        }
        if (i >= NUM_BANK) return false;

        if (first[idx1 + idx2] != 0xffffffff && first[idx3 + idx4] != 0xffffffff) {
            splitting->dst1 = first[idx1 + idx2];
            splitting->dst2 = first[idx3 + idx4];
        }
        else if (first[idx1 + idx3] != 0xffffffff && first[idx2 + idx4] != 0xffffffff) {
            splitting->dst1 = first[idx1 + idx3];
            splitting->dst2 = first[idx2 + idx4];
        }
        else if (first[idx1 + idx4] != 0xffffffff && first[idx2 + idx3] != 0xffffffff) {
            splitting->dst1 = first[idx1 + idx4];
            splitting->dst2 = first[idx2 + idx3];
        }
        else if (first[idx1] != 0xffffffff && first[idx2 + idx3 + idx4] != 0xffffffff) {
            splitting->dst1 = first[idx1];
            splitting->dst2 = first[idx2 + idx3 + idx4];
        }
        else if (first[idx2] != 0xffffffff && first[idx1 + idx3 + idx4] != 0xffffffff) {
            splitting->dst1 = first[idx2];
            splitting->dst2 = first[idx1 + idx3 + idx4];
        }
        else if (first[idx3] != 0xffffffff && first[idx1 + idx2 + idx4] != 0xffffffff) {
            splitting->dst1 = first[idx3];
            splitting->dst2 = first[idx1 + idx2 + idx4];
        }
        else if (first[idx4] != 0xffffffff && first[idx1 + idx2 + idx3] != 0xffffffff) {
            splitting->dst1 = first[idx4];
            splitting->dst2 = first[idx1 + idx2 + idx3];
        }
        else return false;

    }
    splitting->src = segment;
    return true;
}

void merge_Tile(workload *wl, merge *merging, split *splitting, int tidx) {
    start_Time(TIME_MERGE);
    int num_segment = tmat.ptr[tidx + 1] - tmat.ptr[tidx];
    int first_segment = tmat.ptr[tidx];
    int num_merge = 0;
    for (int segment = tmat.ptr[tidx]; segment < tmat.ptr[tidx + 1]; ++segment) {
        uint16_t target = tmat.row_hash[segment];
        if (first[target] == 0xffffffff) {
            first[target] = segment - first_segment;
        }
        else {
            merging[first_segment + num_merge].src = tmat.mapping[segment];
            merging[first_segment + num_merge].dst = first[target];
            tmat.num[segment] = 0;
            num_merge++;
        }
    }

    int cur_segment = num_segment - num_merge;
    int target = cur_segment & 0x0000000f;
    int num_split = 0;
    wl[tidx].num_Split = 0;

    // Comment IF block to test w/o Segment Split
    if (target > 0 && target < THRESHOLD) {
        for (int segment = tmat.ptr[tidx + 1] - 1; segment >= tmat.ptr[tidx]; --segment) {
            if (split_Segments(first, splitting + (tidx * THRESHOLD + num_split), segment, tmat.num[segment])) {
                num_split++;
                if (num_split == target) break;
            }
        }
        if (num_split == target) {
            for (int i = 0; i < num_split; ++i) {
                tmat.num[splitting[tidx * THRESHOLD + i].src] = 0;
                splitted[splitting[tidx * THRESHOLD + i].src] = i;
                splitting[tidx * THRESHOLD + i].src = tmat.mapping[splitting[tidx * THRESHOLD + i].src];
            }
            wl[tidx].num_Split = num_split;
            cur_segment -= wl[tidx].num_Split;
        }
        else num_split = 0;
    }

    int count = 0;
    for (int segment = tmat.ptr[tidx]; segment < tmat.ptr[tidx + 1]; ++segment) {
        if (tmat.num[segment] > 0) {
            mapping[segment] = count;
            tmat.mapping[count] = tmat.mapping[segment];
            tmat.num[count] = tmat.num[segment];
            tmat.row_hash[count] = tmat.row_hash[segment];
            count++;
        }
    }
    for (int merge = 0; merge < num_merge; ++merge) {
        if (splitted[first_segment + merging[first_segment + merge].dst] != 0xff) {
            int index = splitted[first_segment + merging[first_segment + merge].dst];
            merging[first_segment + merge].dst = 0x80000000 | index;
            splitting[tidx * THRESHOLD + index].src = 0x80000000 | index;
        }
        else {
            merging[first_segment + merge].dst = mapping[merging[first_segment + merge].dst];
        }
    }
    for (int split = 0; split < num_split; ++split) {
        splitting[tidx * THRESHOLD + split].dst1 = mapping[splitting[tidx * THRESHOLD + split].dst1];
        splitting[tidx * THRESHOLD + split].dst2 = mapping[splitting[tidx * THRESHOLD + split].dst2];
    }

    wl[tidx].num_Segment = cur_segment;
    wl[tidx].num_GEMM = ((cur_segment + NUM_BANK - 1) >> 4);
    wl[tidx].num_Merge = num_merge;
    for (int segment = tmat.ptr[tidx]; segment < tmat.ptr[tidx + 1]; ++segment) first[tmat.row_hash[segment]] = 0xffffffff;

    end_Time(TIME_MERGE);

    tmat.ne -= (wl[tidx].num_Merge + wl[tidx].num_Split);
}


void merge_Tile_Non_Binary(workload *wl, merge *merging, split *splitting, int tidx) {
    start_Time(TIME_MERGE);
    int num_segment = tmat.ptr[tidx + 1] - tmat.ptr[tidx];
    int first_segment = tmat.ptr[tidx];
    int num_merge = 0;
    for (int segment = tmat.ptr[tidx]; segment < tmat.ptr[tidx + 1]; ++segment) {
        if (tmat.num[segment] < 16) continue;
        uint16_t target = tmat.row_hash[segment];
        if (first[target] == 0xffffffff) {
            first[target] = segment - first_segment;
        }
        else {
            merging[first_segment + num_merge].src = tmat.mapping[segment];
            merging[first_segment + num_merge].dst = first[target];
            tmat.num[segment] = 0;
            num_merge++;
        }
    }

    int cur_segment = num_segment - num_merge;
    int target = cur_segment & 0x0000000f;
    int num_split = 0;
    wl[tidx].num_Split = 0;

    // Comment IF block to test w/o Segment Split
    if (target > 0 && target < THRESHOLD) {
        for (int segment = tmat.ptr[tidx + 1] - 1; segment >= tmat.ptr[tidx]; --segment) {
            if (split_Segments(first, splitting + (tidx * THRESHOLD + num_split), segment, tmat.num[segment])) {
                num_split++;
                if (num_split == target) break;
            }
        }
        if (num_split == target) {
            for (int i = 0; i < num_split; ++i) {
                tmat.num[splitting[tidx * THRESHOLD + i].src] = 0;
                splitted[splitting[tidx * THRESHOLD + i].src] = i;
                splitting[tidx * THRESHOLD + i].src = tmat.mapping[splitting[tidx * THRESHOLD + i].src];
            }
            wl[tidx].num_Split = num_split;
            cur_segment -= wl[tidx].num_Split;
        }
        else num_split = 0;
    }

    int count = 0;
    for (int segment = tmat.ptr[tidx]; segment < tmat.ptr[tidx + 1]; ++segment) {
        if (tmat.num[segment] > 0) {
            mapping[segment] = count;
            tmat.mapping[count] = tmat.mapping[segment];
            tmat.num[count] = tmat.num[segment];
            tmat.row_hash[count] = tmat.row_hash[segment];
            count++;
        }
    }
    for (int merge = 0; merge < num_merge; ++merge) {
        if (splitted[first_segment + merging[first_segment + merge].dst] != 0xff) {
            int index = splitted[first_segment + merging[first_segment + merge].dst];
            merging[first_segment + merge].dst = 0x80000000 | index;
            splitting[tidx * THRESHOLD + index].src = 0x80000000 | index;
        }
        else {
            merging[first_segment + merge].dst = mapping[merging[first_segment + merge].dst];
        }
    }
    for (int split = 0; split < num_split; ++split) {
        splitting[tidx * THRESHOLD + split].dst1 = mapping[splitting[tidx * THRESHOLD + split].dst1];
        splitting[tidx * THRESHOLD + split].dst2 = mapping[splitting[tidx * THRESHOLD + split].dst2];
    }

    wl[tidx].num_Segment = cur_segment;
    wl[tidx].num_GEMM = ((cur_segment + NUM_BANK - 1) >> 4);
    wl[tidx].num_Merge = num_merge;
    for (int segment = tmat.ptr[tidx]; segment < tmat.ptr[tidx + 1]; ++segment) first[tmat.row_hash[segment]] = 0xffffffff;

    end_Time(TIME_MERGE);

    tmat.ne -= (wl[tidx].num_Merge + wl[tidx].num_Split);
}
