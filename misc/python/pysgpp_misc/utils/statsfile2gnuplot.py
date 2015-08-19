# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#!/usr/bin/python
# This file is part of SGpp (pysgpp), a program package making use of spatially adaptive sparse grids to solve numerical problems
# 
# Copyright (C) 2007 Dirk Pflueger (pflueged@in.tum.de)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

## @brief Create gnuplot-plot from statsfile.
#
# Plot e.g. accuracy against lambda and number of grid points
# @version $CURR$
from bin import tools
import optparse, sys

# parse args
parser = optparse.OptionParser()
parser.set_usage('%prog [options] statsfile')
parser.add_option("-o", action="store", type="string", dest="outfile", help="Filename where the calculated alphas are stored")
parser.add_option("-i", action="store", type="int",default="3", metavar="n", dest="ignore", help="Number of first entries with general information before tripels (|grid|, tr_acc, te_acc) start. Default: 3")
parser.add_option("-x", "--x1", action="store", type="int",default="2", metavar="n", dest="x1", help="Number of entry to be used as x1-dimension")
parser.add_option("-t", "--type", action="store", type="choice",default="test", dest="type", choices=["test", "train"], help="Choose whether test or training acc. is considered (def. test)")
(options,args)=parser.parse_args()
if len(args) < 1:
    parser.parse_args(['-h'])

# set offset for test or training acc.
if options.type == "train":
    tt = 1
else:
    tt = 2

# read statsfile
if args[0] == "-" and hasattr(sys.stdin,"readlines"):
    data = sys.stdin.readlines()
else:
    f = open(args[0], 'r')
    data = f.readlines()
    f.close()

## Helper function to sorts data by x1
# @param x the data string
# @return x1 as float
def keyFun(x):
    y = x.strip()
    if y == "":
        return 0
    else:
        return float(y.split(", ")[options.x1-1])
    
data.sort(key=keyFun, reverse=True)

s = '''#set terminal wxt
#set output "bla.png"
set pm3d implicit at s
set hidden3d offset 1 trianglepattern 3 undefined 1 altdiagonal bentover
set style data lines
set contour base
set logscale x
#splot "bupagross.dat" with lines #, "ripley_data.dat" with points
splot "-" with lines
'''

# iterate over data
for line in data:
    line = line.strip()
    if line != "":
        line = line.split(", ")
        if len(line)==1:
            line = line.split(None)
        i = options.ignore
        j = 0
        while i < len(line):
            s += "%s %d %s\n" % (line[options.x1-1],j,line[i+tt])
            i += 3
            j += 1
        s += "\n"

# output
if options.outfile:
    tools.writeStringToFile(s, options.outfile)
else:
    print s