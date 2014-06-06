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
	#  mod linear
	../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} modlinear 4 0.000000316228 250 0.0001 0 0.0 100 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_modlinear_${PREC}_${VEC}_Level_4.log
	../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} modlinear 5 0.00000001 250 0.0001 0 0.0 100 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_modlinear_${PREC}_${VEC}_Level_5.log
	../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} modlinear 6 0.000000316228 250 0.0001 0 0.0 100 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_modlinear_${PREC}_${VEC}_Level_6.log
	../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} modlinear 7 0.000001 250 0.0001 0 0.0 100 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_modlinear_${PREC}_${VEC}_Level_7.log
	../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} modlinear 8 0.000001 250 0.0001 0 0.0 100 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_modlinear_${PREC}_${VEC}_Level_8.log
	../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} modlinear 9 0.000001 250 0.0001 0 0.0 100 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_modlinear_${PREC}_${VEC}_Level_9.log
	../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} modlinear 10 0.000001 250 0.0001 0 0.0 100 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_modlinear_${PREC}_${VEC}_Level_10.log
	../bin/ClassifyBenchmark ../input/chess_05D_small_3fields_tr.dat.arff ../input/chess_05D_small_3fields_te.dat.arff 0 ${PREC} modlinear 8 0.000001 190 0.0001 35 0.0 500 20 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_small_modlinear_${PREC}_${VEC}_adaptive.log
done
