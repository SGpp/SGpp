#!/usr/bin/python

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org


##
# Helper tool, creating projections of grids.
# Creates gnuplot printable projections of a sparse grid provided with --grid.
# Output either as eps, png, or to stdout.
# Useful, e.g., to examine the results of adaptive refinement criteria


import sys, os, re, optparse
sys.path.append(os.environ["SGPP"])#+"/bin")
from bin import tools
from pysgpp import *


if __name__ == '__main__':
    # parse args
    parser = optparse.OptionParser()
    parser.set_usage('''%prog 
         Creates gnuplot printable projections of a sparse grid provided with --grid.
         Output either as eps, png, or to stdout.
    ''')
    parser.add_option("--grid", dest="grid", action="store", type="string", 
                      help="Filename for of grid")
    parser.add_option("--prefix", dest="prefix", action="store", 
                      default=None, 
                      help="Prefix of eps-, png- and gnuplot-files. If not set, no files are created")
    parser.add_option("--eps", dest="eps", action="store_true", 
                      default=False, 
                      help="Create eps")
    parser.add_option("--png", dest="png", action="store_true", 
                      default=False, 
                      help="Create png")
    parser.add_option("--nolabel", dest="nolabel", action="store_true", 
                      default=False, 
                      help="Don't plot label")
    parser.add_option("--monochrome", dest="monochrome", action="store_true", 
                      default=False, 
                      help="Don't use colors. Only for --eps")
    parser.add_option("--blackwhite", dest="blackwhite", action="store_true", 
                      default=False, 
                      help="Plot all grid points the same way")
    parser.add_option("--dotwidth", dest="dotwidth", action="store", type="int",
                      default=None, 
                      help="Width of dots (gnuplot lw)")
    (options,args)=parser.parse_args()
    
    # check arguments
    if not options.grid:
        print "--grid missing"
        parser.parse_args(['-h'])
        
    # read grid
    grid = tools.readGrid(options.grid)
    
    gridStorage = grid.getStorage()
    dim = gridStorage.getDimension()
    p = DataVector(dim)
    
    # prepare dotwidth
    if options.dotwidth:
        dotwidth = " lw %d" % (options.dotwidth)
    else:
        dotwidth = ""
    
    # prepare multiplot
    if True:
        if options.nolabel:
            xoffset = 0.05
            yoffset = 0.05
        else:
            xoffset = 1.0/12.0
            yoffset = 1.0/12.0
        dx = (1.0-xoffset)/dim
        dy = (1.0-yoffset)/dim
        if options.nolabel:
            xoffset = 0
            yoffset = 0
    
        s = """
    # reading grid from %s
    set pm3d map
    set nokey
    unset colorbox
    set size square 1,1
    #set size square %f,%f
    #set rmargin 0 
    #set lmargin 0
    #set tmargin 0 
    #set bmargin 0
    #set multiplot layout %d,%d
    #set size square 2,2
    set multiplot
    #
    """ % (options.grid, 1.0/dim, 1.0/dim, dim, dim)
    
        for i in range(dim):
            for j in range(dim):
                s += "set lmargin 0\nset rmargin 0\nset tmargin 0\nset bmargin 0\n"
                #s += "unset label\nset label \"*(x%d,x%d)*\" at graph 0.6,0.6\n" % (i,j)
                #s += "set xlabel 'x%d'\nset ylabel 'x%d'\n" % (i,j)
                s += "set size %f,%f\n" % (1.4*dx, 1.4*dy)
                s += "set origin %f,%f\n" % (i*dx+xoffset*2/3,j*dy+yoffset*2/3)
                s += "set yrange[%g:%g]\n" % (0,1)
                if i==0 and not options.nolabel:
                    s += "set ytics (%g,%g,%g) format \"%%f\"\n" % (0,0.5,1)
                else:
                    s += "set noytic\n"
                s += "set xrange[%g:%g]\n" % (0,1)
                if j == 0 and not options.nolabel:
                    s += "set xtics (%g,%g,%g) rotate by 90 format \"%%f\" offset 0,-2\n" % (0,0.5,1)
                else:
                    s += "set noxtic\n"
                if dim <= 3:
                    if options.blackwhite:
                        s += "splot '-' with points notitle lc rgbcolor \"black\"%s\n" % (dotwidth)
                    else:
                        s += "splot '-' with points notitle palette%s\n" % (dotwidth)
                else:
                    if options.blackwhite:
                        s += "splot '-' with dots notitle lc rgbcolor \"black\"%s\n" % (dotwidth)
                    else:
                        s += "splot '-' with dots notitle palette%s\n" % (dotwidth)
    
                
                # project grid
                gridpoints = {} # count, how often projected grid point at same place
                for k in xrange(gridStorage.getSize()):
                    gridStorage.getPoint(k).getStandardCoordinates(p)
                    if not gridpoints.has_key((p[i], p[j])):
                        gridpoints[(p[i], p[j])] = 1
                    else:
                        gridpoints[(p[i], p[j])] += 1
                    
                for key in gridpoints.keys():
                    s += "%g %g %g\n" % (key[0], key[1], gridpoints[key])
    
                s += "e\n"
                s += "#\n"
    
        s += "set nomultiplot\n"
    
        if options.eps:
            (cin, couterr) = os.popen2('gnuplot')
            if options.monochrome:
                cin.write("""set terminal postscript monochrome enhanced size 10in,10in font "Arial" 10
    set output '%sprojection.eps'"""%(options.prefix) + s)
            else:
                cin.write("""set terminal postscript color enhanced size 10in,10in font "Arial" 10
    set output '%sprojection.eps'"""%(options.prefix) + s)
            cin.close()
            print couterr.read()
        elif options.png:
            (cin, couterr) = os.popen2('gnuplot')
            cin.write("""set terminal png enhanced small size 1000,1000 enhanced
    set output '%sprojection.png'"""%(options.prefix) + s)
            cin.close()
            print couterr.read()
        else:
            #print "set terminal x11"
            print s
