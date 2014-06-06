#!/bin/sh
# run specific variables, please adjust
export OMP_NUM_THREADS=16
# prefix for result files
RESPRE=SNB-EP_E5-2670
VECTYPE="X86SIMD"
PRECISION="SP DP"

###### DO NOT MODIFY BELOW THIS LINE #######
export KMP_AFFINITY=proclist=[0-15],granularity=thread,explicit,verbose
#export GOMP_CPU_AFFINITY="0-15"
for VEC in $VECTYPE
do
	for PREC in $PRECISION
	do
		# test Checkerboard, linear boundary
		../bin/ClassifyBenchmark ../input/chess_05D_3fields_tr.dat.arff ../input/chess_05D_3fields_te.dat.arff 0 ${PREC} linearboundary 3 0.000001 250 0.0001 6 0.0 100 250 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_${VEC}_linearboundary_${PREC}.log
		# test Checkerboard, mod linear
		../bin/ClassifyBenchmark ../input/chess_05D_3fields_tr.dat.arff ../input/chess_05D_3fields_te.dat.arff 0 ${PREC} modlinear 5 0.0000001 250 0.0001 8 0.0 100 250 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_ChkBrd_05D_${VEC}_modlinear_${PREC}.log
		# test DR5, linear boundary
		../bin/ClassifyBenchmark ../input/DR5_nowarnings_less05_train.arff ../input/DR5_nowarnings_less05_test.arff 1 ${PREC} linearboundary 6 0.00001 200 0.0001 7 0.0 200 150 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_DR5_${VEC}_linearboundary_${PREC}.log
		# test DR5, mod linear
		../bin/ClassifyBenchmark ../input/DR5_nowarnings_less05_train.arff ../input/DR5_nowarnings_less05_test.arff 1 ${PREC} modlinear 6 0.00001 200 0.0001 7 0.0 200 150 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_DR5_${VEC}_modlinear_${PREC}.log
		# test friedman1, linear boundary
		../bin/ClassifyBenchmark ../input/friedman1_90000_train.arff ../input/friedman1_10000_test.arff 1 ${PREC} linearboundary 0 0.0000001 200 0.0001 2 0.0 50 200 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_friedman1_${VEC}_linearboundary_${PREC}.log
		# test friedman1, mod linear
		../bin/ClassifyBenchmark ../input/friedman1_90000_train.arff ../input/friedman1_10000_test.arff 1 ${PREC} modlinear 2 0.000001 200 0.0001 2 0.0 10 200 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_friedman1_${VEC}_modlinear_${PREC}.log
		# test friedman2, linear boundary
		../bin/ClassifyBenchmark ../input/friedman2_90000_train_norm.arff ../input/friedman2_10000_test_norm.arff 1 ${PREC} linearboundary 2 0.00001 200 0.0001 0 0.0 100  200 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_friedman2_${VEC}_linearboundary_${PREC}.log
		# test friedman2, mod linear
		../bin/ClassifyBenchmark ../input/friedman2_90000_train_norm.arff ../input/friedman2_10000_test_norm.arff 1 ${PREC} modlinear 3 0.000001 200 0.0001 0 0.0 100 200 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_friedman2_${VEC}_modlinear_${PREC}.log
		# test friedman3, linear boundary
		../bin/ClassifyBenchmark ../input/friedman3_90000_train_norm.arff ../input/friedman3_10000_test_norm.arff 1 ${PREC} linearboundary 2 0.0000001 200 0.0001 8 0.0 10 200 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_friedman3_${VEC}_linearboundary_${PREC}.log
		# test friedman3, mod linear
		../bin/ClassifyBenchmark ../input/friedman3_90000_train_norm.arff ../input/friedman3_10000_test_norm.arff 1 ${PREC} modlinear 2 0.0001 200 0.0001 7 0.0 10 200 0.0001 ${VEC} 2>&1 | tee ../log/${RESPRE}_friedman3_${VEC}_modlinear_${PREC}.log
	done
done
