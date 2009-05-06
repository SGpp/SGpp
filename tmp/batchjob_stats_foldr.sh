#!/bin/sh

LEVELS="4 5 6 7 8"
LAMBDAS="0.1 0.01 0.001 0.0001 0.000001"
GRIDS="complete-boundary uscaled-boundary no-boundary border"
FOLDLEVEL="10"
TESTDATA="twospirals.wieland.arff.gz"

for level in $LEVELS
do
	for lambda in $LAMBDAS
	do
		for grid in $GRIDS
		do
			STATSFILE="${TESTDATA}.${grid}"
			echo "python classifier.py --data ${TESTDATA} --level $level --zeh laplace --mode foldr --foldlevel $FOLDLEVEL --lambda ${lambda} -i 2000 --${grid}"
			case "$grid" in
  				"no-boundary") python classifier.py --data $TESTDATA --level $level --zeh laplace --mode foldr --foldlevel $FOLDLEVEL --lambda $lambda -i 2000 --stats $STATSFILE;;
  				"uscaled-boundary") python classifier.py --data $TESTDATA --level $level --zeh laplace --mode foldr --foldlevel $FOLDLEVEL --lambda $lambda --$grid -i 2000 --stats $STATSFILE;;
  				"complete-boundary") python classifier.py --data $TESTDATA --level $level --zeh laplace --mode foldr --foldlevel $FOLDLEVEL --lambda $lambda --$grid -i 2000 --stats $STATSFILE;;
  				"border") python classifier.py --data $TESTDATA --level $level --zeh laplace --mode foldr --foldlevel $FOLDLEVEL --lambda $lambda --$grid -i 2000 --stats $STATSFILE;;
  			esac
		done
	done
done