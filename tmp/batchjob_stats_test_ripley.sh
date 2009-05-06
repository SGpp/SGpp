#!/bin/sh

LEVELS="1 2 3 4 5 6 7 8 9"
LAMBDAS="1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.09 0.08 0.07 0.06 0.05 0.04 0.03 0.02 0.01 0.009 0.008 0.007 0.006 0.005 0.004 0.003 0.002 0.001 0.0009 0.0008 0.0007 0.0006 0.0005 0.0004 0.0003 0.0002 0.0001 0.00009 0.00008 0.00007 0.00006 0.00005 0.00004 0.00003 0.00002 0.00001 0.000009 0.000008 0.000007 0.000006 0.000005 0.000004 0.000003 0.000002 0.000001 0.0000009 0.0000008 0.0000007 0.0000006 0.0000005 0.0000004 0.0000003 0.0000002 0.0000001 0.00000001 0.000000001"
GRIDS="complete-boundary uscaled-boundary no-boundary border"
TRAINDATA="../datasets/ripley/ripleyGarcke.train.arff.gz"
TESTDATA="../datasets/ripley/ripleyGarcke.test.arff.gz"
RERS="0.0001 0.00001 0.000001 0.0000001 0.00000001"

for lambda in $LAMBDAS
do
	for level in $LEVELS
	do
		for grid in $GRIDS
		do
			for rer in $RERS
			do
				STATSFILE="${TRAINDATA}.${grid}.level_${level}.r_${rer}"
				echo "python classifier.py --data ${TRAINDATA} --test ${TESTDATA} --level ${level} --zeh laplace --mode test --lambda ${lambda} -i 10000 -r ${rer} --${grid}"
				case "$grid" in
	  				"no-boundary") python classifier.py --data $TRAINDATA --test $TESTDATA --level $level --zeh laplace --mode test --lambda $lambda -i 10000 -r $rer --stats $STATSFILE;;
	  				"uscaled-boundary") python classifier.py --data $TRAINDATA --test $TESTDATA --level $level --zeh laplace --mode test --lambda $lambda --$grid -i 10000 -r $rer --stats $STATSFILE;;
	  				"complete-boundary") python classifier.py --data $TRAINDATA --test $TESTDATA --level $level --zeh laplace --mode test --lambda $lambda --$grid -i 10000 -r $rer --stats $STATSFILE;;
	  				"border") python classifier.py --data $TRAINDATA --test $TESTDATA --level $level --zeh laplace --mode test --lambda $lambda --$grid -i 10000 -r $rer --stats $STATSFILE;;
	  			esac
  			done
		done
	done
done