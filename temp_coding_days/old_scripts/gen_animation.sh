#!/bin/sh

# generate all frames (jpeg) that are needed to build 
# an animation
for file in *.gnuplot
do
	cp $file solution.gnuplot
	gnuplot gen_animation_frame.gplt
	cp temp.jpg ${file}.jpg
	rm solution.gnuplot
	rm temp.jpg
done
# now build an mpeg movie form the jpeg files
ls *.jpg | jpeg2yuv -f 25 -I p | mpeg2enc -o animation.mpg