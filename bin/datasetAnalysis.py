# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

# This file is part of SGClass, a program package making use of spatially adaptive sparse grids to solve numerical problems
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with pyclass. If not, see <http://www.gnu.org/licenses/>.

## @package classifier
# @ingroup bin
# @brief Show some statistics for data files
#Handles either ARFF-files or plain data files. If the ARFF file
#contains a class attribute, it also shows statistics for the class
#distribution. 
#
# Help with <tt>--help</tt>.
# @version $CURR$

from optparse import OptionParser
import sys, os, math
from tools import *

parser = OptionParser()
parser.set_usage('%prog [options]\n\t  Gives some statistics for datasets, given either in arff or simple file format')
parser.add_option("-i", "--infile", action="append", type="string", dest="infiles", help="Specifies the inputfiles to analyse.")
(options,args)=parser.parse_args()

if options.infiles == None:
    parser.parse_args(['-h'])
	

# loop over infiles
for filename in options.infiles:
    try:
        # read in files
        print "================= %20s =================" %(filename)
        ftype = isARFFFile(filename)
        if ftype == ARFF:
            dataset = readDataARFF(filename)
        elif ftype == SIMPLE:
            dataset = readDataTrivial(filename)
        else:
            sys.stderr.write("Skipping "+filename+os.linesep)
            continue

        # analyse data
        # header
        dim = len(dataset["data"])
        numpoints = len(dataset["data"][0])
        print "Dim (#attributes): %d"%(dim)
        print " %-4s %-12s %-12s %-12s %-12s %-12s" % ("Dim",
                                                       "Min",
                                                       "Max",
                                                       "mean",
                                                       "unbiased V",
                                                       "samplestddev")
        # traverse all attributes
        for i in range(dim):
            total_sum = sum(dataset["data"][i])
            mean = total_sum/float(numpoints)
            unbiased_variance = sum(map(lambda x: (x-mean)**2, dataset["data"][i]))/float(numpoints-1)
            sample_stddev = math.sqrt(unbiased_variance)
            print "  %02d %12f %12f %12f %12f %12f" %(i+1,
                                                 min(dataset["data"][i]),
                                                 max(dataset["data"][i]),
                                                 mean,
                                                 unbiased_variance,
                                                 sample_stddev)
        print "#data points: %d"%(numpoints)
        # statistics for class distribution
        if dataset.has_key("classes"):
            print "Class distribution:"
            class_count = {}
            for c in dataset["classes"]:
                if class_count.has_key(c):
                    class_count[c] += 1
                else:
                    class_count[c] = 1
            class_values = class_count.keys()
            class_values.sort()
            for c in class_values:
                print "  %12f %d" % (c, class_count[c])

    except Exception, e:
        sys.stderr.write("ERROR: Skipping "+filename+os.linesep)
        print "  ",e