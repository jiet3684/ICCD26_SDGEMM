#! /bin/bash

SDGEMM_DIR=.
DATASET_DIR=$SDGEMM_DIR/synthetic
# DATASET_DIR=$SDGEMM_DIR/synthetic/skewness
RESULT_FILE=$SDGEMM_DIR/SDGEMM_RESULTS.csv

# TECHNIQUE_LIST=("CPU" "T" "T+M" "T+B" "T+O" "T+M+B" "T+M+O" "T+B+O" "T+M+B+O")
# OPTIMIZATION_LIST=("-B" "-t" "-t -m" "-t -b" "-t -o" "-t -m -b" "-t -m -o" "-t -b -o" "-t -m -b -o")
TECHNIQUE_LIST=("CPU" "T" "T+M" "T+M+B" "T+M+B+O")
OPTIMIZATION_LIST=("-B" "-t" "-t -m" "-t -m -b" "-t -m -b -o")
MATRIX_LIST="4096.10 4096.20 4096.30 4096.40 4096.50 12288.10 12288.20 12288.30 4096.sk1 4096.sk2 4096.sk3 4096.sk4 12288.sk1 12288.sk2 12288.sk3 12288.sk4"
# MATRIX_LIST="4096.sk1 4096.sk2 4096.sk3 4096.sk4 12288.sk1 12288.sk2 12288.sk3 12288.sk4"


rm $RESULT_FILE
touch $RESULT_FILE

# echo ",CPU(FPGA),Uncompressed,Tiled SpMM,Tiled SpMM + M,Tiled SpMM + M + B, Tiled SpMM + M + O, Tiled SpMM + M + B + O" >> $RESULT_FILE

echo ",,Build tCSR,Matrix Compression,Matrix Align,DMA Transfer,GEMM,Overlapped,Total Time,Merge,Split,Segment,GEMM,Orininal Tile,Binded Tile,CSR Size,tCSR Size" >> $RESULT_FILE
for i in $MATRIX_LIST
do
	
    for ((j=0;j<5;j++))
    do
        printf "COMMAND: $SDGEMM_DIR/spmm ${OPTIMIZATION_LIST[$j]} $DATASET_DIR/$i"
        printf "$i,${TECHNIQUE_LIST[$j]}" >> $RESULT_FILE
        $SDGEMM_DIR/spmm -d ${OPTIMIZATION_LIST[$j]} $DATASET_DIR/$i
    done
    echo "" >> $RESULT_FILE
done
