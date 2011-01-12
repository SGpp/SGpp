#!/usr/bin/python
# Copyright (C) 2009 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp

##
# @package data_projections
# @ingroup bin.utils
# @author Dirk Pflueger
# @version $CURR$
# Helper tool, plotting projections of data.
# Plots the grid file provided by --grid

from bin import classifier, tools
import optparse

if __name__ == '__main__':
    # parse args
    parser = optparse.OptionParser()
    parser.set_usage('''%prog 
         Plot function with gnuplot.
    ''')
    parser.add_option("-g", "--gnuplot", action="store", type="string", dest="gnuplot", help="In 2D case, the generated can be stored in a gnuplot readable format.")
    parser.add_option("--grid", action="store", type="string", dest="grid", help="Filename for Grid-resume. For fold? and test. Full filename.")
    parser.add_option("-A", "--alpha", action="store", type="string", dest="alpha", help="Filename for a file containing an alpha-Vector")
    parser.add_option("-R", "--resolution", action="store", type="int",default=50, metavar="RESOLUTION", dest="res", help="Specifies the resolution of the gnuplotfile")
    (options,args)=parser.parse_args()
    
    # check arguments
    if not (options.gnuplot and options.grid and options.alpha):
        print "--file missing"
        parser.parse_args(['-h'])
    
    
    grid = tools.readGrid(options.grid)
    alpha = classifier.buildTrainingVector(classifier.openFile(options.alpha))
    tools.writeGnuplot(options.gnuplot, grid, alpha, options.res)

