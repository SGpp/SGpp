# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

## This is alex's test file

from optparse import OptionParser
import sys
from tools import *
from toolsExtended import *
from pysgpp import *
from painlesscg import cg,sd,cg_new
from math import sqrt
import random
import numpy as np

from array import array

sys.path.append("../bin")

try:
    import psyco
    psyco.full()
    print "Using psyco"
except:
    pass


#-------------------------------------------------------------------------------
## Outputs a deprecated warning for an option
# @param option Parameter set by the OptionParser
# @param opt Parameter set by the OptionParser
# @param value Parameter set by the OptionParser
# @param parser Parameter set by the OptionParser
def callback_deprecated(option, opt, value, parser):
    print "Warning: Option %s is deprecated." % (option)
    
    
#-------------------------------------------------------------------------------
## calculated the condition number of the DM Systemmatrix
def calc_condition():
    factory = Grid.createLinearGrid(6)
    level = 3
    gen = factory.getGenerator()
    gen.regular(level)
    
    training = buildTrainingVector(openFile('datasets/bupa_liver/liver-disorders_normalized.arff.gz'))
    
    aem = 345
    lam = 0.001
    
    print "Number of gridpoints:" + str(factory.getSize())
    print "generating laplacian matrix..."
    laplace_m = generateCMatrix(factory)
    print laplace_m
    print "generating B*B^T matrix..."
    B_res = generateBBTMatrix(factory, training) #np.dot(B_m,Bt_m)
    print B_res
    print "multiplying aem*lambda*C..."
    C = aem * lam * laplace_m
    print "adding C and B_res..."
    C = C + B_res
    print C
    print "calculating condition number..."
    cond = np.linalg.cond(C)
    
    print cond
    
    
#===============================================================================
# Main
#===============================================================================

# check so that file can also be imported in other files
if __name__=='__main__':
    #start the test programm
    calc_condition()