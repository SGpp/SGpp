#!/bin/bash

[[ ! -d "$1" ]] && echo "First argument must be directory." && exit 1

./bin/graph.py "$1/fgrid"
mv graph*.png "$1"

gnuplot -e "file_path='$1/ferr0'" -e "output_file_path='$1/plot1.png'" ./bin/plot
gnuplot -e "file_path1='$1/ferr1'" -e "file_path2='$1/ferr2'" -e "output_file_path='$1/plot2.png'" ./bin/plot2
