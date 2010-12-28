#!/bin/sh

LEVELS="0 1 2 3"
REFINEMENTS="0 1 2 3 4"
LAMBDAS="0.1 0.01 0.001 0.0001 0.00001 0.000001 0.0000001"
REFINEPOINTS="10 50"
CGMAX=200
CGEPS=0.0004
TRAINFILE="friedman1_10000_train.arff"
TESTFILE="friedman1_10000_test.arff"
PREFIXOUT="Friedman1_tr-10000_te-10000"
REFINETHRESH=0.0

for level in $LEVELS
do
	for refNums in $REFINEMENTS
	do
		for lambda in $LAMBDAS
		do
			for numRefPoi in $REFINEPOINTS
			do
			OUTFILE="${PREFIXOUT}_sl-${level}_refCount-${refNums}_lam-${lambda}_refPoints-${numRefPoi}.out"
			echo "./ClassifyBenchmark_OCL_NVOffload $TRAINFILE $TESTFILE 1 $level $lambda $CGMAX $CGEPS $refNums $REFINETHRESH $numRefPoi"
			./ClassifyBenchmark_OCL_NVOffload $TRAINFILE $TESTFILE 1 $level $lambda $CGMAX $CGEPS $refNums $REFINETHRESH $numRefPoi> $OUTFILE
			done
		done
	done
done

