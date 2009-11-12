#!/bin/sh

LEVELS="4 5 6 7 8 9"
LAMBDAS="1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.09 0.08 0.07 0.06 0.05 0.04 0.03 0.02 0.01 0.009 0.008 0.007 0.006 0.005 0.004 0.003 0.002 0.001 0.0009 0.0008 0.0007 0.0006 0.0005 0.0004 0.0003 0.0002 0.0001 0.00005 0.00001 0.000005 0.000001 0.0000005 0.0000001 0.00000001 0.000000001"
GRIDS="complete-boundary trapezoid-boundary no-boundary border"
FOLDLEVEL="10"
TESTDATA="../datasets/twospirals/twospirals.wieland.arff"

for level in $LEVELS
do
	for lambda in $LAMBDAS
	do
		for grid in $GRIDS
		do
			STATSFILE="${TESTDATA}.${level}.${grid}"
			echo "python classifier.py --data ${TESTDATA} --level $level --zeh laplace --mode foldr --foldlevel $FOLDLEVEL --lambda ${lambda} -i 10000 --${grid}"
			case "$grid" in
  				"no-boundary") python classifier.py --data $TESTDATA --level $level --zeh laplace --mode foldr --foldlevel $FOLDLEVEL --lambda $lambda -i 10000 --stats $STATSFILE;;
  				"trapezoid-boundary") python classifier.py --data $TESTDATA --level $level --zeh laplace --mode foldr --foldlevel $FOLDLEVEL --lambda $lambda --$grid -i 10000 --stats $STATSFILE;;
  				"complete-boundary") python classifier.py --data $TESTDATA --level $level --zeh laplace --mode foldr --foldlevel $FOLDLEVEL --lambda $lambda --$grid -i 10000 --stats $STATSFILE;;
  				"border") python classifier.py --data $TESTDATA --level $level --zeh laplace --mode foldr --foldlevel $FOLDLEVEL --lambda $lambda --$grid -i 10000 --stats $STATSFILE;;
  			esac
		done
	done
done