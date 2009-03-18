#!/bin/sh

LEVELS="3 4 5 6 7 8"
LAMBDAS="1 0.1 0.01 0.001 0.0001 0.00001 0.000001"
GRIDS="complete-boundary uscaled-boundary no-boundary"

TESTDATA="twospirals.wieland.arff.gz"

#generate gnuplot data
for level in $LEVELS
do
	for lambda in $LAMBDAS
	do
		for grid in $GRIDS
		do
			GNUFILE="${TESTDATA}.level${level}.lambda${lambda}.${grid}.gnuplot"
			echo "python classifier.py --data ${TESTDATA} --level $level --zeh laplace --mode normal --gnuplot ${GNUFILE} --resolution 200 --lambda ${lambda} -i 2000 --${grid}"
			case "$grid" in
  				"no-boundary") python classifier.py --data $TESTDATA --level $level --zeh laplace --mode normal --gnuplot $GNUFILE --resolution 200 --lambda $lambda -i 2000;;
  				"uscaled-boundary") python classifier.py --data $TESTDATA --level $level --zeh laplace --mode normal --gnuplot $GNUFILE --resolution 200 --lambda $lambda --$grid -i 2000;;
  				"complete-boundary") python classifier.py --data $TESTDATA --level $level --zeh laplace --mode normal --gnuplot $GNUFILE --resolution 200 --lambda $lambda --$grid -i 2000;;
  			esac
		done
	done
done

#generate heatmaps with generated data
for file in *.gnuplot
do
	cp $file temp.gnu
	gnuplot print_gnuplot_heatmaps
	cp temp.png ${file}.png
	rm temp.gnu
	rm temp.png
done
