# This file is part of sgpp, a program package making use of spatially adaptive sparse grids to solve numerical problems.
#
# Copyright (C) 2007  Joerg Blank (blankj@in.tum.de), Richard Roettger (roettger@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with pyclass. If not, see <http://www.gnu.org/licenses/>.
#

## @package classifier
# @ingroup bin
# @brief Main script to do classification, regression, ...
# @version $CURR$


from optparse import OptionParser
import sys
from tools import *
from math import sqrt
import random

from pysgpp import *
from painlesscg import cg,sd,cg_new
from providers.griddecorator import *
from providers.evalproviders import *
from providers.refineproviders import *
from providers.dataproviders import *

try:
    import psyco
    psyco.full()
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
## Formats a list of type mainlist as a string
# <pre>
#     main_list  = {plain_list | string}*
#     plain_list = string*
# </pre>
# @param l Mainlist
def format_optionlist(l):
    def join_inner_list(entry):
        if type(entry) is list:
            return "("+' OR '.join(entry)+")"
        else:
            return entry
    return ' '*4+' AND '.join(map(lambda entry: join_inner_list(entry), l))


#-------------------------------------------------------------------------------
## Checks whether a valid mode is specified,
# whether all required options for the mode are given and
# executes the corresponding action (function)
# @param mode current mode
def exec_mode(mode):

    if mode=="help":
        print "The following modes are available:"
        for m in modes.keys():
            print "%10s: %s" % (m, modes[m]['help'])
        sys.exit(0)
    
    # check valid mode
    if not modes.has_key(mode):
        print("Wrong mode! Please refer to --mode help for further information.")
        sys.exit(1)

    # check if all required options are set
    attributesValid = True
    for requiredAttribute in modes[mode]['required_options']:
        # OR
        if type(requiredAttribute) is list:
            isAtleastOneSet = False
            for alternative in requiredAttribute:
                if getattr(options, alternative, None):
                    isAtleastOneSet = True
            attributesValid = attributesValid and isAtleastOneSet
        else:
            if not getattr(options, requiredAttribute, None):
                attributesValid = False
    if not attributesValid:
        print ("Error!")
        print ("For --mode %s you have to specify the following options:\n" % (mode)
               + format_optionlist(modes[mode]['required_options']))
        sys.exit(1)

    # execute action
    modes[mode]['action'](mode)
    


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def openFile(filename):
    try:
        data = readDataARFF(filename)
    except:
        print ("An error occured while reading " + filename + "!")
        sys.exit(1)
        
    if data.has_key("classes") == False:
        print ("No classes found in the given File " + filename + "!")
        sys.exit(1)
        
    return data

## Load data from file
def loadData():
    data = []
    try:
        for filename in options.data:
            data.append(readDataARFF(filename))
    except:
        raise Exception("An error occured while reading " + filename + "!")

    
    if len(data) < 1:
        raise Exception("No data files supplied.")
    
    return data




#-------------------------------------------------------------------------------
def saveCheckpoint(gridDecorator, checkpoint, adaptStep):
    import pickle
    filename = "%s_%08d" % (checkpoint, gridDecorator.getSize())
    fout = open(filename + ".a" + str(options.adapt_start + adaptStep) + ".grid", "wb")
    pickle.dump(gridDecorator, fout, -1)
    fout.close()

#===============================================================================
# General algorithms
#===============================================================================


def learningStep(gridDecorator, data):
    if gridDecorator.eval:
        gridDecorator.eval.reset()
    
    if gridDecorator.refine:
        gridDecorator.refine.reset()
    #TODO: Exception, what if eval or refine are None?
        
    for trainingData, testingData in DataProvider(data, gridDecorator.mode, options):
        alpha = DataVector(gridDecorator.getSize())
        alpha.setAll(0.0)

        m = Matrix(gridDecorator.grid, trainingData[0], options.l, gridDecorator.zeh, gridDecorator.basetype)
        b = m.generateb(trainingData[1])

        res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, False, options.verbose)
        print res
        
        if gridDecorator.eval:
            gridDecorator.eval.updateResults(alpha, trainingData, testingData)
        
        if gridDecorator.refine:
            gridDecorator.refine.add(alpha, trainingData)
        
        #writeAlphaARFF(filename + ".alpha.arff", alpha)


def learningAlgorithm(mode):
    data = loadData()
    (data, dim) = DataProvider(data, mode, options).construct()

    gridDecorator = GridDecorator.create(options, dim, mode)
    gridDecorator.setRefine(options.refine)
    gridDecorator.setEval(options.eval)
    
    learningStep(gridDecorator, data)
    if gridDecorator.eval:
        gridDecorator.eval.updateStatistics(options)
    
    #perform adaptive steps if necessary
    for adaptStep in xrange(options.adaptive):
        
        if gridDecorator.refine:
            gridDecorator.refine.refine(options)
            
        print "adaptStep: %d" %adaptStep
        learningStep(gridDecorator, data)
        
        if options.checkpoint:
            saveCheckpoint(gridDecorator, options.checkpoint, adaptStep)
        
        if gridDecorator.eval:
            gridDecorator.eval.updateStatistics(options)
    

def evalAlgorithm(mode):
    
    raise Exception("eval unsupported.")
    
    data = loadData()
    data, dim = DataProvider(data, mode, options).construct()

    gridDecorator = GridDecorator.create(options, dim, mode)

    alpha = DataVector(gridDecorator.getSize())
    #Fill alpha Vector
    
    for training, testing in data_providers[gridDecorator.mode]["split"](data):
        gridDecorator.eval.testing(alpha, testing)
        
    gridDecorator.eval.updateStatistics(options)
    


#===============================================================================
# Main
#===============================================================================

#-------------------------------------------------------------------------------
def doUnsupported(mode):
    raise NotImplementedError(mode + " mode unsupported")

# specifiy the modes:
# modes is an array containing all modes, the options needed by the mode and the action
# that is to be executed
modes = {
    'apply'  : {'help': "classify a dataset with an existing grid (compute accuracy)",
                'required_options': ['data', 'alpha', ['level', 'grid']],
                'action': doUnsupported},
    'eval'   : {'help': "evaluate a sparse grid function given by grid and alphas at data points",
                'required_options': ['data', 'alpha', ['level', 'grid']],
                'action': doUnsupported},
    'normal' : {'help': "learn a dataset",
                'required_options': ['data', ['level', 'grid']],
                'action': learningAlgorithm},
    'test'   : {'help': "learn a dataset with a test dataset",
                'required_options': ['data', 'test'],
                'action': doUnsupported},
    'fold'   : {'help': "learn a dataset with a random n-fold",
                'required_options': ['data'],
                'action': learningAlgorithm},
    'folds'  : {'help': "learn a dataset with a sequential n-fold",
                'required_options': ['data'],
                'action': learningAlgorithm},
    'foldr'  : {'help': "learn a dataset with a stratified n-fold",
                'required_options': ['data'],
                'action': learningAlgorithm},
    'foldf'  : {'help': "learn a dataset with a n-fold from a set of files",
                'required_options': ['data'],
                'action': doUnsupported}
    }
def _main():
    # Initialize OptionParser, set Options
    parser = OptionParser()
    parser.add_option("-l", "--level", action="store", type="int", dest="level", help="Gridlevel")
    parser.add_option("-a", "--adaptive", action="store", type="int", default="0", dest="adaptive", metavar="NUM", help="Using an adaptive Grid with NUM of refines")
    parser.add_option("--adapt_points", action="store", type="int", default="1", dest="adapt_points", metavar="NUM", help="Number of points in one refinement iteration")
    parser.add_option("--adapt_start", action="store", type="int", default="0", dest="adapt_start", metavar="NUM", help="The index of adapt step to begin with")
    parser.add_option("-m", "--mode", action="store", type="string", default="apply", dest="mode", help="Specifies the action to do. Get help for the mode please type --mode help.")
    parser.add_option("-C", "--zeh", action="store", type="string", default="laplace", dest="zeh", help="Specifies the action to do.")
    parser.add_option("-f", "--foldlevel", action="store", type="int",default="10", metavar="LEVEL", dest="f_level", help="If a fold mode is selected, this specifies the number of sets generated")
    parser.add_option("-L", "--lambda", action="store", type="float",default="0.000001", metavar="LAMBDA", dest="l", help="Lambda")
    parser.add_option("-i", "--imax", action="store", type="int",default="400", metavar="MAX", dest="imax", help="Max number of iterations")
    parser.add_option("-r", "--accuracy", action="store", type="float",default="0.0001", metavar="ACCURACY", dest="r", help="Specifies the accuracy of the CG-Iteration")
    parser.add_option("-d", "--data", action="append", type="string", dest="data", help="Filename for the Datafile.")
    parser.add_option("-t", "--test", action="store", type="string", dest="test", help="File containing the testdata")
    parser.add_option("-A", "--alpha", action="store", type="string", dest="alpha", help="Filename for a file containing an alpha-Vector")
    
    parser.add_option("-b", "--base", action="store", type="string", default="linear", dest="basetype", help="Specifies basetype. Values: linear, modlinear, poly, modpoly")
    parser.add_option("-e", "--eval", action="store", type="string", default="classes", dest="eval", help="Specifies evaluation method. Values: classes, regression")
    parser.add_option("-R", "--refine", action="store", type="string", default="surplus", dest="refine", help="Specifies refinement method. Values: surplus, error")
    
    parser.add_option("-P", "--polynom", action="store", type="int", default="-1", dest="polynom", help="Specifies maximal polynomial degree for polynomial base functions")
    
    parser.add_option("-s", "--stats", action="store", type="string", dest="stats", help="In this file the statistics from the test are stored")
    parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose", help="Provides extra output")

    parser.add_option("--seed", action="store", type="int", dest="seed", help="Random seed used for initializing")
    parser.add_option("--checkpoint", action="store", type="string", dest="checkpoint", help="Filename for checkpointing. For fold? and test. No file extension.")
    parser.add_option("--grid", action="store", type="string", dest="grid", help="Filename for Grid-resume. For fold? and test. Full filename.")
    # parse options
    global options
    options,args = parser.parse_args()
    
    # check some options
    zeh = options.zeh.lower()
    
    options.basetype = options.basetype.lower()
    options.eval = options.eval.lower()
    options.refine = options.refine.lower()
    
    # check C-mode
    if zeh == "help":
        print "The following C-modes are available:"
        for m in zeh_modes.keys():
            print "%15s: %s" % (m, zeh_modes[m])
        sys.exit(0)
    elif zeh not in zeh_modes.keys():
        print("Wrong C-mode! Please refer to -C help for further information.")
        sys.exit(1)

    # Execute the mode
    exec_mode(options.mode.lower())



if __name__=="__main__":
    _main()
