#!/bin/gnuplot

#set terminal pngcairo size 800,600 enhanced font "Verdana,10"
set mouse
set xrange [0:1]
set yrange [0:1]
set zrange [-2:2]
set cbrange [-1:1] 
set xlabel "x"
set ylabel "y"
set style line 1 lt -1 lw 0.3
set pm3d hidden3d 1
do for [i = 1:100] {
    splot "out/solution.dat" index (i-1) using 1:2:3 with pm3d
    pause 0.02
}
reread
pause mouse keypress
