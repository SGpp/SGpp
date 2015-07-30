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


## @package statsfileInfo
# @ingroup bin
# @brief Compute statistics for statsfile
# @version $CURR$

import optparse, math, sys, tools

# parse args
parser = optparse.OptionParser()
parser.set_usage('%prog [options] statsfile\n   Compute statistics for statsfile')
parser.add_option("-o", action="store", type="string", dest="outfile", help="Write it to a file rather than stdout")
parser.add_option("-i", action="store", type="int",default="3", metavar="n", dest="ignore", help="Number of first entries with general information before tripels (|grid|, tr_acc, te_acc) start. Default: 3")
parser.add_option("-x", "--x1", action="store", type="int",default="2", metavar="n", dest="x1", help="Number of column with lambda. Default: 2")
parser.add_option("--lambdas", action="store_true", default=True, dest="lambdas", help="Show statistics for every lambda.")
parser.add_option("--summary", action="store_true", default=True, dest="summary", help="Show summary")
parser.add_option("--min", action="store_true", default=False, dest="min", help="Consider minimum rather than maximum. Default: maximum.")
(options,args)=parser.parse_args()
if len(args) < 1:
    parser.parse_args(['-h'])

# read statsfile
if args[0] == "-" and hasattr(sys.stdin,"readlines"):
    txt = sys.stdin.readlines()
else:
    f = open(args[0], 'r')
    txt = f.readlines()
    f.close()

# parse data: ommit empty lines, split and convert to float
data = []
for i in range(len(txt)):
    txt[i] = txt[i].strip()
    if txt[i] != "":
        data.append(map(lambda x: float(x), txt[i].split(", ")))

# sort data by lambda
def sortFun(a,b):
    if a[options.x1-1] < b[options.x1-1]:
        return -1
    elif a[options.x1-1] == b[options.x1-1]:
        return 0
    else:
        return 1
data.sort(sortFun, reverse=True)

# init values
if options.min:
    # minimum
    opt_tr = sys.maxint
    opt_te = sys.maxint
else:
    # maximum
    opt_tr = 0
    opt_te = 0
opt_tr_numpoints = sys.maxint
opt_te_numpoints = sys.maxint
opt_tr_s = ""
opt_te_s = ""
s = ""

if options.lambdas:
    s = "Lambda      counter - |grid|   Train     Test     - |grid|   Train     Test\n"

# iterate over data
for line in data:
        # get min and max
        if options.min:
            # minimum
            lopt_tr = sys.maxint
            lopt_te = sys.maxint
        else:
            # maximum
            lopt_tr = 0
            lopt_te = 0
        lopt_tr_s = ""
        lopt_te_s = ""
        i = options.ignore
        while i < len(line):
            # training
            if ((not options.min and line[i+1] > lopt_tr) or (options.min and line[i+1] < lopt_tr)):
                lopt_tr = line[i+1]
                lopt_tr_numpoints = int(line[i])
                lopt_tr_s = "%7d, %-8g, %-8g" % (
                    int(line[i]), line[i+1], line[i+2])
            # testing
            if ((not options.min and line[i+2] > lopt_te) or (options.min and line[i+2] < lopt_te)):
                lopt_te = line[i+2]
                lopt_te_numpoints = int(line[i])
                lopt_te_s = "%7d, %-8g, %-8g" % (
                    int(line[i]), line[i+1], line[i+2])
            i += 3
        # output
        lambda_s = "%-11g" %(line[options.x1-1])
        counter = "%06.2f" %(math.log(line[options.x1-1], 0.1**(1/30.0)))
        if options.lambdas:
            s = (s + lambda_s + "  " + counter + " - " + lopt_tr_s + " - " + lopt_te_s + "\n")
        if ((not options.min and lopt_tr > opt_tr) or (options.min and lopt_tr < opt_tr) or (lopt_tr == opt_tr and lopt_tr_numpoints < opt_tr_numpoints)):
            opt_tr = lopt_tr
            opt_tr_numpoints = lopt_tr_numpoints
            opt_tr_s = lambda_s + ", " + lopt_tr_s
        if ((not options.min and lopt_te > opt_te) or (options.min and lopt_te < opt_te) or (lopt_te == opt_te and lopt_te_numpoints < opt_te_numpoints)):
            opt_te = lopt_te
            opt_te_numpoints = lopt_te_numpoints
            opt_te_s = lambda_s + ", " + lopt_te_s
            
if options.summary:
    if options.lambdas:
        s = s + "\n"
    s = (s + "Lambda      |grid|   Train     Test\n" +
         opt_tr_s + "\n"+ opt_te_s)
            
# output
s = "Results for file %s\n" % (args[0]) + s

if options.outfile:
    tools.writeStringToFile(s, options.outfile)
else:
    print s
