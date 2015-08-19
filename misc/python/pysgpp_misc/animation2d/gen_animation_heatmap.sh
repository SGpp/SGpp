#!/bin/sh

# generate all frames (jpeg) that are needed to build 
# an animation
for file in *.gnuplot
do
	cp $file solution.gnuplot
	gnuplot gen_animation_frame_heatmap.gplt
	cp temp.png ${file}.png
	rm solution.gnuplot
	rm temp.png
done
# now build an mpeg movie form the png files
mogrify *.png
png2yuv -j %014d.gnuplot.png -f25 -Ip -L0 | mpeg2enc -f12 -a3 -x1280 -y720 -o animation.mpg
