#! /bin/bash

SDGEMM_DIR=.
DATASET_DIR=$SDGEMM_DIR/../dataset
RESULT_FILE=$SDGEMM_DIR/SDGEMM_RESULTS.csv

TECHNIQUE_LIST=("T" "T+M" "T+M+B" "T+M+B+O")
OPTIMIZATION_LIST=("-t" "-t -m" "-t -m -b" "-t -m -b -o")
#MATRIX_LIST="nd3k.mtx gupta3.mtx TSOPF_RS_b2383_c1.mtx thread.mtx bcsstk32.mtx mycielskian16.mtx crankseg_1.mtx gupta2.mtx kron_g500-logn16.mtx nd24k.mtx consph.mtx TSOPF_FS_b300_c3.mtx"
MATRIX_LIST="crankseg_1.mtx"


rm $RESULT_FILE
touch $RESULT_FILE

# echo ",CPU(FPGA),Uncompressed,Tiled SpMM,Tiled SpMM + M,Tiled SpMM + M + B, Tiled SpMM + M + O, Tiled SpMM + M + B + O" >> $RESULT_FILE

echo ",,Build tCSR,Matrix Compression,Matrix Align,DMA Transfer,GEMM,Overlapped,Total Time,Merge,Split,Segment,GEMM,Orininal Tile,Binded Tile,CSR Size,tCSR Size" >> $RESULT_FILE
for i in $MATRIX_LIST
do
    for ((j=0;j<4;j++))
    # for ((j=0;j<9;j++))
    do
        printf "COMMAND: $SDGEMM_DIR/spmm ${OPTIMIZATION_LIST[$j]} $DATASET_DIR/$i"
        printf "$i,${TECHNIQUE_LIST[$j]}" >> $RESULT_FILE
        $SDGEMM_DIR/spmm -d ${OPTIMIZATION_LIST[$j]} $DATASET_DIR/$i
    done
    echo "" >> $RESULT_FILE
done
