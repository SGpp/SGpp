#!/bin/sh
# run specific variables, please adjust
export OMP_NUM_THREADS=24
# prefix for result files
RESPRE=WSM-EP
VECTYPE="SSE"
PRECISION="SP DP"
LEVELS="2 3 4 5 6"
REFINEMENTS="0 1 2 3 4 5 6 7"

# best known config, so far: Level 5, 5 refinements

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
				../bin/ClassifyBenchmark ../input/cod-rna_train_unit.arff ../input/cod-rna_test_unit.arff 0 ${PREC} linearboundary ${LEVEL} 0.000001 200 0.0001 ${REFINE} 0.0 100 40 0.01 2>&1 | tee ../log/${RESPRE}_cod-rna_${VEC}_linearboundary_${PREC}_Level_${LEVEL}_Refinements_${REFINE}.log
				#  mod linear
				../bin/ClassifyBenchmark ../input/cod-rna_train_unit.arff ../input/cod-rna_test_unit.arff 0 ${PREC} modlinear ${LEVEL} 0.000001 200 0.0001 ${REFINE} 0.0 100 40 0.01 2>&1 | tee ../log/${RESPRE}_cod-rna_${VEC}_modlinear_${PREC}_Level_${LEVEL}_Refinements_${REFINE}.log
			done
		done
	done
done
