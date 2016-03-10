#!/usr/bin/python

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org


import sys, os, re, optparse, random, math
sys.path.append("../../bin/")
import tools
from pysgpp import *

random.seed(12345)
N = 10000

X = DataMatrix(N, 11)
p = DataVector(11)

for i in xrange(N):
    for d in range(10):
        if d == 0 or d == 1:
            p[d] = random.uniform(-0.5,0.5)
        else:
            p[d] = random.random()

    eps = random.normalvariate(0.0,0.1)
    p[10] = p[0]**2 * p[1]**2 + eps

    p[0] += 0.5
    p[1] += 0.5

    X.setRow(i, p)

fd = tools.gzOpen("parabola.10000_test.arff", 'w')

fd.write("""@RELATION "Parabola %s"\n\n""" % (N))
for d in range(X.getNcols()-1):
    fd.write("""@ATTRIBUTE x%d NUMERIC\n""" % (d))
fd.write("""@ATTRIBUTE class NUMERIC\n\n@DATA\n""")
for i in xrange(X.getNrows()):
    X.getRow(i, p)
    fd.write(','.join([str(p[d]) for d in range(X.getNcols())])+"\n")

fd.close()