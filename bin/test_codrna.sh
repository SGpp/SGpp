#!/bin/sh
# run specific variables, please adjust
export OMP_NUM_THREADS=24
# prefix for result files
RESPRE=WSM-EP
VECTYPE="SSE"
PRECISION="SP DP"
LEVELS="2 3 4"
REFINEMENTS="0 1 2 3 4 5 6 7"

###### DO NOT MODIFY BELOW THIS LINE #######
export KMP_AFFINITY=compact,granularity=thread,verbose
. /opt/intel/bin/compilervars.sh intel64

for VEC in $VECTYPE
do
	for PREC in $PRECISION
	do
		for LEVEL in $LEVELS
		do
			for REFINE in $REFINEMENTS
			do
				#  linear boundary
				#../bin/ClassifyBenchmark ../input/CR_train_unit.labeled ../input/CR_test_unit.labeled 0 ${PREC} linearboundary ${LEVEL} 0.000001 40 0.01 ${REFINE} 0.0 100 200 0.0001 2>&1 | tee ../log/${RESPRE}_cod-rna_${VEC}_linearboundary_${PREC}_Level_${LEVEL}_Refinements_${REFINE}.log
				#  mod linear
				../bin/ClassifyBenchmark ../input/CR_train_unit.labeled ../input/CR_test_unit.labeled 0 ${PREC} modlinear ${LEVEL} 0.000001 40 0.01 ${REFINE} 0.0 100 200 0.0001 2>&1 | tee ../log/${RESPRE}_cod-rna_${VEC}_modlinear_${PREC}_Level_${LEVEL}_Refinements_${REFINE}.log
			done
		done
	done
done
