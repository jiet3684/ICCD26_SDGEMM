#include "SpMM.hpp"

csr mat;
tcsr tmat;
void free_Memory(int mode);
int parse_option(int argc, char * const argv[]);
// #define POWER

int main(int argc, char **argv) {
    int mode;
    if ((mode = parse_option(argc, argv)) == -1) return 0;

    // convert_CSR(argv[argc - 1]);
	convert_synthetic(argv[argc - 1]);


    if (check_option(mode, BASE_SPMM)) {
        SpMM_CPU();
        GEMM_CPU();
        print_Result(mode);
        free_Memory(mode);
        return 0;
    }


    uint32_t A_ROW_PAD = tmat.nr * REG_SIZE; // M, p
    uint32_t A_COL_PAD = (mat.nc + REG_SIZE - 1) / REG_SIZE * REG_SIZE; // K, q

    VAL_TYPE *B = (VAL_TYPE*)malloc(A_COL_PAD * COL_B * sizeof(VAL_TYPE));

    VAL_TYPE *align_A;
    VAL_TYPE *align_B;
    if (check_option(mode, OVERLAPPING)) { // For Double Buffering
        align_A = (VAL_TYPE*)malloc(NUM_BANK * A_COL_PAD * sizeof(VAL_TYPE) << 1);
        align_B = (VAL_TYPE*)malloc(A_COL_PAD * COL_B * sizeof(VAL_TYPE) << 1);
    }
    else {
        align_A = (VAL_TYPE*)malloc(NUM_BANK * A_COL_PAD * sizeof(VAL_TYPE));
        align_B = (VAL_TYPE*)malloc(A_COL_PAD * COL_B * sizeof(VAL_TYPE));
    }
    VAL_TYPE *C = (VAL_TYPE*)malloc(A_ROW_PAD * COL_B * sizeof(VAL_TYPE));
    
#ifdef POWER
    long cpu = 0, dram = 0;
    long temp_uj;
    FILE *fp;
    fp = fopen("/sys/class/powercap/intel-rapl/intel-rapl:0/energy_uj", "rt");
    fscanf(fp, "%ld", &temp_uj);
    cpu -= temp_uj;
    fclose(fp);
    fp = fopen("/sys/class/powercap/intel-rapl/intel-rapl:1/energy_uj", "rt");
    fscanf(fp, "%ld", &temp_uj);
    cpu -= temp_uj;
    fclose(fp);
    fp = fopen("/sys/class/powercap/intel-rapl/intel-rapl:0/intel-rapl:0:0/energy_uj", "rt");
    fscanf(fp, "%ld", &temp_uj);
    dram -= temp_uj;
    fclose(fp);
    fp = fopen("/sys/class/powercap/intel-rapl/intel-rapl:1/intel-rapl:1:0/energy_uj", "rt");
    fscanf(fp, "%ld", &temp_uj);
    dram -= temp_uj;
    fclose(fp);
#endif
    
    for (int i = 0; i < ITER; ++i) {
        make_TiledCSR(mode);

        tiled_spmm(align_A, B, align_B, C, A_ROW_PAD, A_COL_PAD, (uint32_t)COL_B, mode);
    }

#ifdef POWER
    fp = fopen("/sys/class/powercap/intel-rapl/intel-rapl:0/energy_uj", "rt");
    fscanf(fp, "%ld", &temp_uj);
    cpu += temp_uj;
    fclose(fp);
    fp = fopen("/sys/class/powercap/intel-rapl/intel-rapl:1/energy_uj", "rt");
    fscanf(fp, "%ld", &temp_uj);
    cpu += temp_uj;
    fclose(fp);
    fp = fopen("/sys/class/powercap/intel-rapl/intel-rapl:0/intel-rapl:0:0/energy_uj", "rt");
    fscanf(fp, "%ld", &temp_uj);
    dram += temp_uj;
    fclose(fp);
    fp = fopen("/sys/class/powercap/intel-rapl/intel-rapl:1/intel-rapl:1:0/energy_uj", "rt");
    fscanf(fp, "%ld", &temp_uj);
    dram += temp_uj;
    fclose(fp);

    long total_uj = cpu;
    double total_time = print_Result(mode) * 1e3;
    double cpu_avg_watt = (((cpu - dram) / total_time) - 55.601664) / ITER;
    double dram_avg_watt = ((dram / total_time) - 13.296205) / ITER;
    double total_avg_watt = cpu_avg_watt + dram_avg_watt;
    printf("Time:\n%lf\nGFLOPS:\n%lf\n", total_time / 1000.0, 2.0 * mat.ne * COL_B / total_time * 1e-3);
    printf("Total Average Watt: \n%lf\n\n", total_avg_watt);
    printf("Total DRAM Watt: \n%lf\n\n", dram_avg_watt);
    printf("Total CPU Watt: \n%lf\n\n", cpu_avg_watt);
    fp = fopen("profile.csv", "at");
    fprintf(fp, "%lf,%lf,%lf\n", 2.0 * mat.ne * COL_B / total_time * 1e-3, dram_avg_watt, cpu_avg_watt);
    fclose(fp);
#else
    double total_time = print_Result(mode);
    FILE *fp = fopen("profile.csv", "at");
    fprintf(fp, "%lf,%lf\n", total_time, 2.0 * mat.ne * COL_B / total_time * 1e-6);
    fclose(fp);
#endif

    free(B);
    free(align_A);
    free(align_B);
    free(C);
    free_Memory(mode);
}

int parse_option(int argc, char * const argv[]) {
    if (argc == 2) {
        puts("Usage: ./spmm (options) input_graph\n\t-B: Baseline SpMM\n\t-u: Uncompressed GEMM\n\t-c: CSR Format\n\t-t: tCSR Format\n\t-m: Segment Merging\n\t-b: Tile Binding\n\t-o: CPU-PIM Overlapping\n\t-p: Parallel Pre-processing\n");
        return -1;
    }
    int mode = 0;
    int c;
    while ((c = getopt(argc, argv, "Buctmopdbpah")) != -1) {
        switch (c)
        {
            case 'B':
                mode = mode | BASE_SPMM;
                break;
            case 'u':
                mode = mode | UNCOMPRESSED;
                break;
            case 'c':
                mode = mode | CSR_SPMM;
                break;
            case 't':
                mode = mode | TCSR_SPMM;
                break;
            case 'm':
                mode = mode | MERGING;
                break;
            case 'o':
                mode = mode | OVERLAPPING;
                break;
            case 'b':
                mode = mode | BINDING;
                break;
            case 'p':
                mode = mode | PARALLEL;
                break;
            case 'a':
                mode = mode |TCSR_SPMM | MERGING | OVERLAPPING | BINDING | PARALLEL;
                break;
            case 'd':
                mode = mode | DEBUGGING;
                break;
            case 'h':
            default:
                puts("Usage: ./spmm (options) input_graph\n\t-B: Baseline SpMM\n\t-u: Uncompressed GEMM\n\t-c: CSR Format\n\t-t: tCSR Format\n\t-m: Segment Merging\n\t-b: Tile Binding\n\t-o: CPU-PIM Overlapping\n");
                return -1;
        }
    }
    return mode;
}

void free_Memory(int mode) {
    free(mat.nnz);
    free(mat.ptr);
    free(mat.idx);
    if (check_option(mode, TCSR_SPMM)) {
        free(tmat.ptr);
        free(tmat.mapping);
        free(tmat.num);
        free(tmat.row_hash);
    }
    if (check_option(mode, CSR_SPMM)) {
        free(tmat.ptr);
        // free(tmat.mapping);
    }
} 
