#! /bin/bash

SDGEMM_DIR=.
# DATASET_DIR=$SDGEMM_DIR/../dataset
# MATRIX_LIST="nd3k.mtx gupta3.mtx TSOPF_RS_b2383_c1.mtx thread.mtx bcsstk32.mtx mycielskian16.mtx crankseg_1.mtx gupta2.mtx kron_g500-logn16.mtx nd24k.mtx consph.mtx TSOPF_FS_b300_c3.mtx"

DATASET_DIR=$SDGEMM_DIR/../synthetic
MATRIX_LIST="4096.sp1 4096.sp2 4096.sp3 4096.sp4 4096.sp5 12288.sp1 12288.sp2 12288.sp3 12288.sp4 12288.sp5 4096.sk1 4096.sk2 4096.sk3 4096.sk4 4096.sk5 12288.sk1 12288.sk2 12288.sk3 12288.sk4 12288.sk5"


RESULT_FILE=$SDGEMM_DIR/profile.csv
rm $RESULT_FILE
touch $RESULT_FILE


for i in $MATRIX_LIST
do
    ./spmm -a $DATASET_DIR/$i.unique
    # ./spmm -t -m -o -b $DATASET_DIR/$i.unique
done
