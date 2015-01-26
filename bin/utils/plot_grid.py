#!/usr/bin/python

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#
# Helper tool, plotting projections of data.
# Plots the grid file provided by --grid

import sys, os
# adapt path so that bin can be found
sys.path.append(os.path.join(sys.path[0], "../.."))
from bin import classifier, tools
import optparse

if __name__ == '__main__':
    # parse args
    parser = optparse.OptionParser()
    parser.set_usage('''%prog 
Evaluate a sparse grid function, provided by grid and surplus vector
on a regular, full grid. The results are written to a file that
can be used by gnuplot.

Example usage for gnuplot:
  1d: plot "filename" with lines
  2d: splot "filename" with pm3d
    ''')
    parser.add_option("-g", "--gnuplot", action="store", type="string", dest="gnuplot", metavar="FILENAME",
                      help="In 2D case, the generated can be stored in a gnuplot readable format.")
    parser.add_option("--grid", action="store", type="string", dest="grid", metavar="FILENAME",
                      help="Filename for Grid-resume. For fold? and test. Full filename.")
    parser.add_option("-A", "--alpha", action="store", type="string", dest="alpha", metavar="FILENAME",
                      help="Filename for a file containing an alpha-Vector")
    parser.add_option("-R", "--resolution", action="store", type="int",default=50, metavar="N", dest="res", help="Specifies the sampling resolution (def. 50)")
    (options,args)=parser.parse_args()
    
    # check arguments
    if not (options.gnuplot and options.grid and options.alpha):
        print "param missing"
        parser.parse_args(['-h'])
    
    
    grid = tools.readGrid(options.grid)
    alpha = tools.readAlpha(options.alpha)
    tools.writeGnuplot(options.gnuplot, grid, alpha, options.res)
