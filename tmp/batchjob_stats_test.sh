#!/bin/sh

LEVELS="2 3 4 5 6 7 8"
LAMBDAS="1 0.1 0.01 0.001 0.0001 0.000001"
GRIDS="complete-boundary uscaled-boundary no-boundary"
TRAINDATA="ripleyGarcke.train.arff.gz"
TESTDATA="ripleyGarcke.test.arff.gz"

for level in $LEVELS
do
	for lambda in $LAMBDAS
	do
		for grid in $GRIDS
		do
			STATSFILE="${TRAINDATA}.${grid}"
			echo "python classifier.py --data ${TRAINDATA} --test ${TESTDATA} --level $level --zeh laplace --mode test --lambda ${lambda} -i 2000 --${grid}"
			case "$grid" in
  				"no-boundary") python classifier.py --data $TRAINDATA --test $TESTDATA --level $level --zeh laplace --mode test --lambda $lambda -i 2000 --stats $STATSFILE;;
  				"uscaled-boundary") python classifier.py --data $TRAINDATA --test $TESTDATA --level $level --zeh laplace --mode test --lambda $lambda --$grid -i 2000 --stats $STATSFILE;;
  				"complete-boundary") python classifier.py --data $TRAINDATA --test $TESTDATA --level $level --zeh laplace --mode test --lambda $lambda --$grid -i 2000 --stats $STATSFILE;;
  			esac
		done
	done
done