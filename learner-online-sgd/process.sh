#!/bin/bash

[[ ! -d "$1" ]] && echo "First argument must be directory." && exit 1

./bin/grid.py "$1/fgrid"
mkdir -p "$1/grids"
mv *.png "$1/grids"

gnuplot -e "file_path0='$1/ferr0'" -e "file_path1='$1/ferr1'" -e "file_path2='$1/ferr2'" -e "output_file_path='$1/plot-error.png'" ./bin/plot-error
