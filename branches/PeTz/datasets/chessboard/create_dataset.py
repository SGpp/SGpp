#!/usr/bin/python
# -*- coding: latin-1 -*

#############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
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
