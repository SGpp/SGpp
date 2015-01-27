#!/bin/sh

#generate heatmaps with generated data
for file in *.gnuplot
do
	cp $file temp.gnu
	gnuplot print_gnuplot_heatmaps
	cp temp.png ${file}.png
	rm temp.gnu
	rm temp.png
done