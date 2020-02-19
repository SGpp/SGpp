#!/usr/bin/python

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import random, math

numElements = 500000
fileName = "friedman2_4d_" + str(numElements) + ".arff"
seed = 135

# fixed
dim = 4

random.seed(seed)

dat = ""
for i in xrange(numElements):
    p = []
    # dim 1
    p.append(random.uniform(0, 100))
    # dim 2
    p.append(random.uniform(40*math.pi, 560*math.pi))
    # dim 2
    p.append(random.uniform(0, 1))
    # dim 2
    p.append(random.uniform(1, 11))
    eps = random.normalvariate(0.0,125.0)
#     eps = 0.0
    # $ \left( x_0^2 + \left( x_1 x_2 - (x_1 x_3)^{-1} \right)^2 \right)^0.5 + \epsilon $
    # target value
    p.append((p[0]**2 + (p[1]*p[2] - 1.0/(p[1]*p[3]))**2 )**0.5 + eps)
    
    # normalize dataset to [0,1]^d
    p[0] = (p[0] - 0) / (100 - 0)
    p[1] = (p[1] - (40 * math.pi)) / ((560*math.pi) - (40 * math.pi))
    p[2] = (p[2] - 0) / (1 - 0)
    p[3] = (p[3] - 1) / (11 - 1)
    
    for d in xrange(dim + 1):
        if d > 0:
            dat += ","
        dat += str(p[d])
    dat += "\n"
    
outFile = open(fileName, "w")

header = "@RELATION \"" + fileName + "\"\n"
header += "\n"
for d in range(dim):
    header += "@ATTRIBUTE x" + str(d) + " NUMERIC\n"
header += "@ATTRIBUTE class NUMERIC"
header += "\n"
header += "@DATA\n"

outFile.write(header + dat)

print "all done!"