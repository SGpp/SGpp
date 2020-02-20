#!/usr/bin/python

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org


# double f(std::vector<double> point) {
#     return 16.0 * (point[0] - 1) * point[0] * (point[1] - 1) * point[1];
# }

import random
import operator

samples = 100000
dim = 4
fileName = "parabola_simple_" + str(dim) + "d.arff"

def parabola(point):
    return reduce(operator.mul, [-4.0 * (point[d] - 1.0) * point[d] for d in range(len(point))])

dat = ""
for i in range(samples):
    point = []
    for d in range(dim):
        point.append(random.uniform(0.0, 1.0))
        if d > 0:
            dat += ", "
        dat += str(point[d])
        
    dat += ", " + str(parabola(point)) + "\n"

outFile = open(fileName, "w")

header = "@RELATION \"" + fileName + "\"\n"
header += "\n"
for d in range(dim):
    header += "@ATTRIBUTE x" + str(d) + " NUMERIC\n"
header += "@ATTRIBUTE class NUMERIC"
header += "\n"
header += "@DATA\n"

outFile.write(header + dat)
outFile.close()