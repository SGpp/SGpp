#!/bin/sh
# run specific variables, please adjust
export OMP_NUM_THREADS=24
# prefix for result files
RESPRE=WSM-EP
VEC="X86SIMD"
PRECISION="SP DP"

###### DO NOT MODIFY BELOW THIS LINE #######
export KMP_AFFINITY=compact,granularity=thread,verbose
. /opt/intel/bin/compilervars.sh intel64
for PREC in $PRECISION
do
	# linear boundary
	../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} linearboundary 4 0.000000316228 250 0.0001 0 0.0 100 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_linearboundary_${PREC}_${VEC}_Level_4.log
	../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} linearboundary 5 0.00000001 250 0.0001 0 0.0 100 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_linearboundary_${PREC}_${VEC}_Level_5.log
	../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} linearboundary 6 0.000000316228 250 0.0001 0 0.0 100 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_linearboundary_${PREC}_${VEC}_Level_6.log
	../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} linearboundary 7 0.000001 250 0.0001 0 0.0 100 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_linearboundary_${PREC}_${VEC}_Level_7.log
	#../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} linearboundary 8 0.000001 250 0.0001 0 0.0 100 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_linearboundary_${PREC}_${VEC}_Level_8.log
	#../bin/ClassifyBenchmark .../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} linearboundary 9 0.000001 250 0.0001 0 0.0 100 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_linearboundary_${PREC}_${VEC}_Level_9.log
	../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} linearboundary 5 0.000001 200 0.0001 14 0.0 475 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_linearboundary_${PREC}_${VEC}_adaptive.log
done
