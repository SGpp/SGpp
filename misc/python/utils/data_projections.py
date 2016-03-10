#!/usr/bin/python

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org


#
# Helper tool, creating projections of data.
# Creates gnuplot projections of a data set provided with --data.
# Output either as eps, png, or to stdout.
#

import sys, os, re, optparse
if os.environ.has_key("SGPP"):
    sys.path.append(os.environ["SGPP"]+"/bin")
from bin import tools


def getminmax(filename, separator=None):
    """Computes min and max for each parameter"""
    fd = tools.gzOpen(filename, 'r')

    mmax = [-10e10 for d in range(dim)]
    mmin = [10e10 for d in range(dim)]
    for line in fd.readlines():
    	if "@" in line or len(line)==0:
    		continue
        line = line.split(separator)
        for d in range(len(line)-1):
            mmax[d] = max(float(line[d]), mmax[d])
            mmin[d] = min(float(line[d]), mmin[d])
    fd.close()

    return (mmin, mmax)



if __name__=='__main__':
    # parse args
    parser = optparse.OptionParser()
    parser.set_usage('''%prog 
         Creates gnuplot printable projections of a sparse grid provided with --grid.
         Output either as eps, png, or to stdout.
    ''')
    parser.add_option("--data", dest="data", action="store", type="string", 
                      help="Filename for of data")
    parser.add_option("--datsep", dest="datsep", action="store", type="string", 
                      help="Data file separator (default: ' ')", default=" ")
    parser.add_option("--minima", dest="minima", action="store", type="string", 
                      help="Minima for intervals as a string with one entry per attribute, e.g. \"0.1 0.0 0.3 0.3\" in four dimensions", 
                      default=None)
    parser.add_option("--maxima", dest="maxima", action="store", type="string", 
                      help="Maxima for intervals as a string with one entry per attribute, e.g. \"0.9 1.0 0.8 0.7\" in four dimensions", 
                      default=None)
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
    parser.add_option("--dotwidth", dest="dotwidth", action="store", type="int",
                      default=None, 
                      help="Width of dots (gnuplot lw)")
    (options,args)=parser.parse_args()
    
    # check arguments
    if not options.data:
        print "--data missing"
        parser.parse_args(['-h'])
        
    # read data
    fd = tools.gzOpen(options.data, 'r')
    data = fd.readlines()
    fd.close()
    
    if not options.datsep.strip():
        options.datsep = None
        gnuplot_separator = ""
    else:
        gnuplot_separator = 'set datafile separator "%s"' % (options.datsep)
    
    
    if options.data[-3:] == ".gz":
        fname = "< zcat "+options.data
    else:
        fname = options.data
    # get dimension:
    # khakhutv: since at the beginning of the document there can be header, it makes
    # sense to take one the last lines
    line = data[-2].strip().split(options.datsep)
    dim = len(line)-1
    print "# dim =", dim
    
    # prepare dotwidth
    if options.dotwidth:
        dotwidth = " lw %d" % (options.dotwidth)
    else:
        dotwidth = ""
    
    # get bounds of intervals
    if options.minima and options.maxima:
        mmin = map(lambda x: float(x), options.minima.split())
        mmax = map(lambda x: float(x), options.maxima.split())
    else:
        (mmin, mmax) = getminmax(options.data, options.datsep)
    
    # prepare multiplot
    if True:
        margin = 0 #0.05
        if options.nolabel:
            xoffset = 0.02
            yoffset = 0.02
        else:
            xoffset = 1.0/12.0
            yoffset = 1.0/12.0
        dx = (1.0-xoffset)/dim
        dy = (1.0-yoffset)/dim
    
        s = """
    # taking data from %s
    %s
    set pm3d map
    set nokey
    unset colorbox
    set size 1,1 #square 
    set rmargin 0 
    set lmargin 0
    set tmargin 0 
    set bmargin 0
    #set multiplot layout 3,3
    set multiplot
    #
    """ % (options.data, gnuplot_separator)
    
        for i in range(dim):
            for j in range(dim):
                s += "set lmargin 0\nset rmargin 0\nset tmargin 0\nset bmargin 0\n"
                s += "set size %f,%f\n" % (1.38*dx-margin, 1.38*dy-margin)
                s += "set origin %f,%f\n" % (i*dx+xoffset-0.16*dx,j*dy+yoffset-0.16*dy)
                s += "set yrange[%g:%g]\n" % (mmin[j], mmax[j]) # Enforce certain max and min values
                if i==0 and not options.nolabel:
                    s += "set ytics nomirror #(%g,%g,%g) format \"%%11.4e\"\n" % (0,0.5,1)
                else:
                    s += "set noytic\n"
                s += "set xrange[%g:%g]\n" % (mmin[i], mmax[i]) # Enforce certain max and min values
                if j == 0 and not options.nolabel:
                    s += "set xtics nomirror rotate #(%g,%g,%g) format \"%%f\" offset 0,-2\n" % (0,0.5,1)
                else:
                    s += "set noxtic\n"
                if dim <= 3:
                    s += "splot '%s' using %d:%d:%d with points notitle palette%s\n" % (fname, i+1, j+1, dim+1, dotwidth)
                else:
                    s += "splot '%s' using %d:%d:%d with dots notitle palette%s\n" % (fname, i+1, j+1, dim+1, dotwidth)
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
    #        print "set terminal x11"
            print s
