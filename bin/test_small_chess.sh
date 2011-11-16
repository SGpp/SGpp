#!/bin/sh
# run specific variables, please adjust
export OMP_NUM_THREADS=8
# prefix for result files
RESPRE=SNB
VECTYPE="AVX"
PRECISION="SP DP"

###### DO NOT MODIFY BELOW THIS LINE #######
export KMP_AFFINITY=compact,granularity=thread,verbose
. /opt/intel/bin/compilervars.sh intel64
for PREC in $PRECISION
do
	#  mod linear
	../bin/ClassifyBenchmark ../input/CB5d_train_unit.labeled ../input/CB5d_test_unit.labeled 0 ${PREC} modlinear 4 0.000000316228 200 0.0001 0 0.0 100 200 0.0001 2>&1 | tee ../log/${RESPRE}_chess_small_modlinear_${PREC}_Level_4.log
	../bin/ClassifyBenchmark ../input/CB5d_train_unit.labeled ../input/CB5d_test_unit.labeled 0 ${PREC} modlinear 5 0.00000001 200 0.0001 0 0.0 100 200 0.0001 2>&1 | tee ../log/${RESPRE}_chess_small_modlinear_${PREC}_Level_5.log
	../bin/ClassifyBenchmark ../input/CB5d_train_unit.labeled ../input/CB5d_test_unit.labeled 0 ${PREC} modlinear 6 0.000000316228 200 0.0001 0 0.0 100 200 0.0001 2>&1 | tee ../log/${RESPRE}_chess_small_modlinear_${PREC}_Level_6.log
	../bin/ClassifyBenchmark ../input/CB5d_train_unit.labeled ../input/CB5d_test_unit.labeled 0 ${PREC} modlinear 7 0.000001 200 0.0001 0 0.0 100 200 0.0001 2>&1 | tee ../log/${RESPRE}_chess_small_modlinear_${PREC}_Level_7.log
	../bin/ClassifyBenchmark ../input/CB5d_train_unit.labeled ../input/CB5d_test_unit.labeled 0 ${PREC} modlinear 8 0.000001 200 0.0001 0 0.0 100 200 0.0001 2>&1 | tee ../log/${RESPRE}_chess_small_modlinear_${PREC}_Level_8.log
done
