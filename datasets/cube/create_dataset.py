#!/usr/bin/python
# -*- coding: latin-1 -*-
import sys, re, os, optparse, random, tools

# parse args
parser = optparse.OptionParser()
parser.set_usage('''%prog [options]\n\tCreate a d-dimensional dataset, uniformly distributed points. All points with class -1 are in a subcube with volume 1/2, originating at the origin. ''')
parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="Output filename. If not specified, output to stdout")
parser.add_option("-d", "--dimension", action="store", type="int", dest="dim", default="2", help="dimensionality of feature space (number of attributes)")
parser.add_option("-n", "--num", action="store", type="int", dest="num", default="10", help="number of grid points")
parser.add_option("--seed", action="store", type="int", dest="seed", default="732402928", help="Random seed used for initializing")
(options,args)=parser.parse_args()
#if len(args) == 0:
#    parser.parse_args(['-h'])

random.seed(options.seed)

# compute threshold
threshold = 0.5**(1.0/float(options.dim))
print "threshold = %f"%(threshold)

dat = ""
for n in range(options.num):
    point = ""
    b = 0
    for i in range(options.dim):
        r = random.random()
        if r >= threshold:
            b = 1
        point = point + "%f "%(r)
        
    if b:
        point = point + "+1\n"
    else:
        point = point + "-1\n"
    dat = dat + point

# output
if options.outfile:
    tools.writeStringToFile(dat, options.outfile)
else:
    print dat
