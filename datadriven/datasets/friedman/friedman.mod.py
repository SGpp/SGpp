# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

#!/usr/bin/python

import sys, os, re, optparse, random, math
sys.path.append("../../bin/")
import tools

from pysgpp import *

# parse args
parser = optparse.OptionParser()
parser.set_usage('''%prog options
     Creates data for friedman1-3
''')
parser.add_option("-f", "--friedman", dest="friedman", action="store", type="int",
                  default=None, 
                  help="Which friedman dataset to create (1-3)")
parser.add_option("--seed", dest="seed", action="store", type="int",
                  default=None, 
                  help="Seed for random generator (optional)")
parser.add_option("-o", "--outfile", dest="outfile", action="store",
                  default=None, 
                  help="Filename of the outfile. Otherwise output to stdout.")
parser.add_option("-N", dest="N", action="store", type="int",
                  default=None, 
                  help="Number of data points to create")
parser.add_option("--type", dest="t", action="store", type="str",default="", help="For Friedman 1. Possible values: UNIFORM, BALL, CUBE")

(options,args)=parser.parse_args()

# check arguments
if not options.N:
    print "-N missing"
    parser.parse_args(['-h'])
if not options.friedman or options.friedman < 1 or options.friedman > 3:
    print "-f missing or wrong"
    parser.parse_args(['-h'])
    

# init
if options.seed:
    random.seed(options.seed)

# friedman1:
if options.friedman == 1:
    namestring = 'Friedman1, %d data points' % (options.N)
    X = DataMatrix(options.N, 11)
    p = DataVector(11)
    dds = range(10)

    for i in xrange(options.N):

        for d in dds:
            p[d] = random.random()

        # Uniform

        if options.t == "UNIFORM":
            pass

        if options.t == "BALL":

            #
            # (p[0], ..., p[4]) are the coordinates of a point in the ball
            #
            r = 0.3

            p[0] = random.uniform(-r, r)
            r_1 = math.sqrt(r ** 2 - p[0] ** 2)

            p[1] = random.uniform(-r_1, r_1)
            r_2 = math.sqrt(r ** 2 - p[0] ** 2 - p[1] ** 2)

            p[2] = random.uniform(-r_2, r_2)
            r_3 = math.sqrt(r ** 2 - p[0] ** 2 - p[1] ** 2 - p[2] ** 2)

            p[3] = random.uniform(-r_3, r_3)
            r_4 = math.sqrt(r ** 2 - p[0] ** 2 - p[1] ** 2 - p[2] ** 2 - p[3] ** 2)

            p[4] = r_4

            if random.choice([1,0]) == 1:
                p[4] *= (-1)

            for d in range(5):
               p[d] += 0.5
        
        if options.t == "CUBE":
            r = 0.2
            axis = random.randrange(0,5)

            for d in range(5):
                if d == axis:
                    p[d] = random.choice([1,-1]) * r
                else:
                    p[d] = random.uniform(-r, r)

                p[d] += 0.5

        eps = random.normalvariate(0.0,1.0)
        # $ 10\sin(\pi x_0 x_1) + 20(x_2-0.5)^2 + 10x_3 + 5x_4 + \epsilon $

        p[10] = 10.0*math.sin(math.pi*p[0]*p[1]) + 20.0*(p[2]-0.5)**2 + 10*p[3] + 5*p[4] + eps
        X.setRow(i,p)

# friedman2:
elif options.friedman == 2:
    namestring = 'Friedman2, %d data points' % (options.N)
    X = DataMatrix(options.N, 5)
    p = DataVector(5)
    for i in xrange(options.N):
        p[0] = random.uniform(0, 100)
        p[1] = random.uniform(40*math.pi, 560*math.pi)
        p[2] = random.uniform(0, 1)
        p[3] = random.uniform(1, 11)
        eps = random.normalvariate(0.0,125.0)
        # $ \left( x_0^2 + \left( x_1 x_2 - (x_1 x_3)^{-1} \right)^2 \right)^0.5 + \epsilon $
        p[4] = ( p[0]**2 + (p[1]*p[2] - 1.0/(p[1]*p[3]))**2 )**0.5 + eps
        X.setRow(i,p)
# friedman3:
elif options.friedman == 3:
    namestring = 'Friedman3, %d data points' % (options.N)
    X = DataMatrix(options.N, 5)
    p = DataVector(5)
    for i in xrange(options.N):
        p[0] = random.uniform(0, 100)
        p[1] = random.uniform(40*math.pi, 560*math.pi)
        p[2] = random.uniform(0, 1)
        p[3] = random.uniform(1, 11)
        eps = random.normalvariate(0.0,0.1)
        # $ \atan \left( \left( x_1 x_2 - (x_1 x_3)^{-1} \right) / x_0 \right) + \epsilon $
        p[4] = math.atan( (p[1]*p[2] - 1.0/(p[1]*p[3])) / p[0] ) + eps
        X.setRow(i,p)
else:
    sys.exit(1)
 
if options.outfile and ".csv" in options.outfile:
    from numpy import savetxt
    #header = ','.join(['x%d'%i for i in xrange(X.getNcols()-1)] + ['classes'])
    savetxt(options.outfile, X.array(), fmt='%.12f', delimiter = ',') 
    sys.exit(1)      
elif options.outfile:
    fd = tools.gzOpen(options.outfile, 'w')
else:
    fd = sys.stdout

fd.write("""@RELATION "%s"\n\n""" % (namestring))
for d in range(X.getNcols()-1):
    fd.write("""@ATTRIBUTE x%d NUMERIC\n""" % (d))
fd.write("""@ATTRIBUTE class NUMERIC\n\n@Data\n""")
for i in xrange(X.getNrows()):
    X.getRow(i, p)
    fd.write(','.join([str(p[d]) for d in range(X.getNcols())])+"\n")
if options.outfile:
    fd.close()
