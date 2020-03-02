#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import random

numElements = 500000
dim = 4
fields = 5
seed = 135

fileName = "chess_" + str(dim) + "d_" + str(numElements) + ".arff"

random.seed(seed)

# compute intervalsize
intervalsize = 1.0/fields
print "intervalsize = %f"%(intervalsize)

dat = ""
for n in range(numElements):
    point = ""
    b = 1.0
    for i in range(dim):
        r = random.random()
        
        rang = r/intervalsize
        irang = int(rang+1.0)
        
        # calculate class
        if (irang % 2) == 0:
            b *= 1.0
        else:
            b *= -1.0
             
        point = point + str(r) + ", "
        
    point = point + str(b) + "\n"
    dat = dat + point

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
