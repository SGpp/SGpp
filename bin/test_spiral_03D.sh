#!/bin/sh
# run specific variables, please adjust
export OMP_NUM_THREADS=24
# prefix for result files
RESPRE=WSM-EP
VECTYPE="X86SIMD"
PRECISION="SP DP"

###### DO NOT MODIFY BELOW THIS LINE #######
export KMP_AFFINITY=compact,granularity=thread,verbose
. /opt/intel/bin/compilervars.sh intel64
for PREC in $PRECISION
do
	#  mod linear
	../bin/ClassifyBenchmark ../input/spiral_03D_train_unit.arff ../input/spiral_03D_test_unit.arff 0 ${PREC} modlinear 4 0.000001 250 0.00000001 0 0.0 100 20 0.000001 ${VEC} 2>&1 | tee ../log/${RESPRE}_spiral_03D_modlinear_${PREC}_Level_4.log
	../bin/ClassifyBenchmark ../input/spiral_03D_train_unit.arff ../input/spiral_03D_test_unit.arff 0 ${PREC} modlinear 5 0.000001 250 0.00000001 0 0.0 100 20 0.000001 ${VEC} 2>&1 | tee ../log/${RESPRE}_spiral_03D_modlinear_${PREC}_Level_5.log
	../bin/ClassifyBenchmark ../input/spiral_03D_train_unit.arff ../input/spiral_03D_test_unit.arff 0 ${PREC} modlinear 6 0.000001 250 0.00000001 0 0.0 100 20 0.000001 ${VEC} 2>&1 | tee ../log/${RESPRE}_spiral_03D_modlinear_${PREC}_Level_6.log
	../bin/ClassifyBenchmark ../input/spiral_03D_train_unit.arff ../input/spiral_03D_test_unit.arff 0 ${PREC} modlinear 7 0.000001 250 0.00000001 0 0.0 100 20 0.000001 ${VEC} 2>&1 | tee ../log/${RESPRE}_spiral_03D_modlinear_${PREC}_Level_7.log
	../bin/ClassifyBenchmark ../input/spiral_03D_train_unit.arff ../input/spiral_03D_test_unit.arff 0 ${PREC} modlinear 8 0.000001 250 0.00000001 0 0.0 100 20 0.000001 ${VEC} 2>&1 | tee ../log/${RESPRE}_spiral_03D_modlinear_${PREC}_Level_8.log
	../bin/ClassifyBenchmark ../input/spiral_03D_train_unit.arff ../input/spiral_03D_test_unit.arff 0 ${PREC} modlinear 9 0.000001 250 0.00000001 0 0.0 100 20 0.000001 ${VEC} 2>&1 | tee ../log/${RESPRE}_spiral_03D_modlinear_${PREC}_Level_9.log
	../bin/ClassifyBenchmark ../input/spiral_03D_train_unit.arff ../input/spiral_03D_test_unit.arff 0 ${PREC} modlinear 10 0.000001 250 0.00000001 0 0.0 100 20 0.000001 ${VEC} 2>&1 | tee ../log/${RESPRE}_spiral_03D_modlinear_${PREC}_Level_10.log
	../bin/ClassifyBenchmark ../input/spiral_03D_train_unit.arff ../input/spiral_03D_test_unit.arff 0 ${PREC} modlinear 4 0.000001 250 0.00000001 7 0.0 100 20 0.000001 ${VEC} 2>&1 | tee ../log/${RESPRE}_spiral_03D_modlinear_${PREC}_adaptive.log
done
