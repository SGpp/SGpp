#!/usr/bin/python

import random

fileName = "chess_3d.arff"
numElements = 100000
dim = 3
fields = dim
seed = 135

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
        
        if i < 2:
            rang = r/intervalsize
            irang = int(rang+1.0)
        
            # calcualte class
            if (irang % 2) == 0:
                b *= 1.0
            else:
                b *= (-1.0)
             
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