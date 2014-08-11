#!/bin/bash


[[ $# -eq 0 ]] && echo "Need at least one directory." && exit 1

for dir
do
  [[ ! -d "$dir" ]] && echo "Argument must be directory." && exit 1
  rm -rf "$dir/grids"
  mkdir -p "$dir/grids"
  ./bin/grid.py "$dir/fgrid" "$dir/grids"
  gnuplot -e "file_path0='$dir/ferr0'" -e "file_path1='$dir/ferr1'" -e "file_path2='$dir/ferr2'" -e "output_file_path='$dir/plot-error.png'" ./bin/plot-error
done

