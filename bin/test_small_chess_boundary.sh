#!/bin/sh
# run specific variables, please adjust
export OMP_NUM_THREADS=24
# prefix for result files
RESPRE=WSM-EP
VECTYPE="SSE"
PRECISION="DP"

###### DO NOT MODIFY BELOW THIS LINE #######
export KMP_AFFINITY=compact,granularity=thread,verbose
. /opt/intel/bin/compilervars.sh intel64
for PREC in $PRECISION
do
	#  mod linear
	../bin/ClassifyBenchmark ../input/CB5d_train_unit.labeled ../input/CB5d_test_unit.labeled 0 ${PREC} linearboundary 4 0.000000316228 250 0.00000001 0 0.0 100 20 0.000001 2>&1 | tee ../log/${RESPRE}_chess_small_linearboundary_${PREC}_Level_4.log
	../bin/ClassifyBenchmark ../input/CB5d_train_unit.labeled ../input/CB5d_test_unit.labeled 0 ${PREC} linearboundary 5 0.00000001 250 0.00000001 0 0.0 100 20 0.000001 2>&1 | tee ../log/${RESPRE}_chess_small_linearboundary_${PREC}_Level_5.log
	../bin/ClassifyBenchmark ../input/CB5d_train_unit.labeled ../input/CB5d_test_unit.labeled 0 ${PREC} linearboundary 6 0.000000316228 250 0.00000001 0 0.0 100 20 0.000001 2>&1 | tee ../log/${RESPRE}_chess_small_linearboundary_${PREC}_Level_6.log
	../bin/ClassifyBenchmark ../input/CB5d_train_unit.labeled ../input/CB5d_test_unit.labeled 0 ${PREC} linearboundary 7 0.000001 250 0.00000001 0 0.0 100 20 0.000001 2>&1 | tee ../log/${RESPRE}_chess_small_linearboundary_${PREC}_Level_7.log
	../bin/ClassifyBenchmark ../input/CB5d_train_unit.labeled ../input/CB5d_test_unit.labeled 0 ${PREC} linearboundary 8 0.000001 250 0.00000001 0 0.0 100 20 0.000001 2>&1 | tee ../log/${RESPRE}_chess_small_linearboundary_${PREC}_Level_8.log
	../bin/ClassifyBenchmark ../input/CB5d_train_unit.labeled ../input/CB5d_test_unit.labeled 0 ${PREC} linearboundary 9 0.000001 250 0.00000001 0 0.0 100 20 0.000001 2>&1 | tee ../log/${RESPRE}_chess_small_linearboundary_${PREC}_Level_9.log
	../bin/ClassifyBenchmark ../input/CB5d_train_unit.labeled ../input/CB5d_test_unit.labeled 0 ${PREC} linearboundary 4 0.000001 250 0.00000001 7 0.0 100 20 0.000001 2>&1 | tee ../log/${RESPRE}_chess_small_linearboundary_${PREC}_adaptive.log
done
