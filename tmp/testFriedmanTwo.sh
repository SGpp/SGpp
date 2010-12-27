#!/bin/sh

LEVELS="2 3 4 5"
REFINEMENTS="0 1 2 3 4 5 6 7"
LAMBDAS="0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.00005 0.00001 0.000005 0.000001 0.0000005 0.0000001"
REFINEPOINTS="10 50 100 150"
CGMAX=200
CGEPS=0.0004
TRAINFILE="friedman2_10000_train_norm.arff"
TESTFILE="friedman2_10000_test_norm.arff"
PREFIXOUT="Friedman2_tr-10000_te-10000"
REFINETHRESH=0.0

#Do tests without adaptivity during solving with classic adativity
for level in $LEVELS
do
	for refNums in $REFINEMENTS
	do
		for lambda in $LAMBDAS
		do
			for numRefPoi in $REFINEPOINTS
			do
			OUTFILE="${PREFIXOUT}_sl-${level}_refCount-${refNums}_lam-${lambda}_refPoints-${numRefPoi}.out"
			echo "./ClassifyBenchmark_OCL_NVOffload $TRRAINFILE $TESTFILE 1 $level $lambda $CGMAX $CGEPS $refNums $REFINETHRESH $numRefPoi"
			./ClassifyBenchmark_OCL_NVOffload $TRRAINFILE $TESTFILE 1 $level $lambda $CGMAX $CGEPS $refNums $REFINETHRESH $numRefPoi> $OUTFILE
			done
		done
	done
done

