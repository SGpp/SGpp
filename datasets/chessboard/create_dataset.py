# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#!/usr/bin/python
# -*- coding: latin-1 -*

#############################################################################
                                    #
#############################################################################

import sys, re, os, optparse, random
sys.path.append("../../bin")
import tools

# parse args
parser = optparse.OptionParser()
parser.set_usage('''%prog [options]\n\tCreate a d-dimensional dataset, uniformly distributed points that looks like a chessboard. ''')
parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="Output filename. If not specified, output to stdout")
parser.add_option("-d", "--dimension", action="store", type="int", dest="dim", default="2", help="dimensionality of feature space (number of attributes)")
parser.add_option("-n", "--num", action="store", type="int", dest="num", default="10", help="number of grid points")
parser.add_option("-f", "--fields", action="store", type="float", dest="fields", default="4.0", help="Number of chess field in one dimension")
parser.add_option("--seed", action="store", type="int", dest="seed", default="732402928", help="Random seed used for initializing")
(options,args)=parser.parse_args()

random.seed(options.seed)

# compute intervalsize
intervalsize = 1.0/options.fields
print "intervalsize = %f"%(intervalsize)

dat = ""
for n in range(options.num):
    point = ""
    b = 1.0
    for i in range(options.dim):
        r = random.random()
        rang = r/intervalsize
        irang = int(rang+1.0)
        
        # calcualte class
        if (irang % 2) == 0:
            b *= 1.0
        else:
            b *= (-1.0)
             
        point = point + str(r) + " "
        
    point = point + str(b) + "\n"
    
    dat = dat + point

# output
if options.outfile:
    tools.writeStringToFile(dat, options.outfile)
else:
    print dat