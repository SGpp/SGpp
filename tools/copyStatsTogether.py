# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#############################################################################
                                    #
#############################################################################

import sys
from math import sqrt
import random
import re
from array import array

def copyStatsTogether():
    f = open('ripleyGarcke.train.arff.gz.no-boundary.level_1.r_0.0001.stats')
    line = f.readlines()
    f.close()
    lines = len(line)
    
    numlevel = 9
    levels = [1,2,3,4,5,6,7,8,9]
    fileTogether = "ripleyGarcke.train.arff.gz.no-boundary.level_All.r_0.0001.stats"
    fileInPre = "ripleyGarcke.train.arff.gz.no-boundary.level_"
    fileInPost = ".r_0.0001.stats"
    fout = file(fileTogether, "w")
    
    for i in range(lines):
        for j in range(numlevel):
            filein = fileInPre + str(levels[j]) + fileInPost
            fin = open(filein, "r")
            for n in range(i):
                fin.readline()
                
            newEntry = fin.readline().strip()
            fin.close()
            fout.write(newEntry)
            fout.write(" ")
            
        fout.write("\n")
    
    
    fout.close()

if __name__=='__main__':
    copyStatsTogether()