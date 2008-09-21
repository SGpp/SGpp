# This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems.
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
from pysgpp import *
from painlesscg import cg,sd,cg_new
from math import sqrt

import random

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
    a = True
    for attrib in modes[mode]['required_options']:
        # OR
        if type(attrib) is list:
            b = False
            for attrib2 in attrib:
                if getattr(options, attrib2, None):
                    b = True
            a = a and b
        else:
            if not getattr(options, attrib, None):
                a = False
    if not a:
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


#-------------------------------------------------------------------------------
## Constructs a new grid.
# If options.grid is set, then read in a stored grid. If not, construct a new
# grid dependent on the dimension dim, on options.level and options.polynom.
# Sets the use of boundary functions according to options.border.
# @param dim the grid dimension
# @return a grid
def constructGrid(dim):
    factories = {"linear" : lambda : Grid.createLinearGrid(dim),
             "modlinear" : lambda : Grid.createModLinearGrid(dim),
             "poly" : lambda : Grid.createPolyGrid(dim, options.polynom),
             "modpoly" : lambda : Grid.createModPolyGrid(dim, options.polynom),
             }
    grid = factories[options.basetype]()
    gen = grid.createGridGenerator()
    gen.regular(options.level)
    
    return grid

#-------------------------------------------------------------------------------
## Construction from serialized grid or new
def construction(dim = None, mode = None):
    if options.grid:
        import pickle
        fin = open(options.grid, "rb")
        
        status = pickle.load(fin)
        
        fin.close()
        return status
    else:
        return Status(constructGrid(dim), mode, options.basetype, options.zeh)

#-------------------------------------------------------------------------------
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
# For n-fold-cv:
# assemble training vector for fold <omit>
#-------------------------------------------------------------------------------
def assembleTrainingVector(dvecs,cvecs,omit):
    size = 0
    for dataset in dvecs:
        size = size + dataset.getSize()
    
    size = size - dvecs[omit].getSize()
    
    training = DataVector(size, dvecs[0].getDim())
    classes = DataVector(size)
    
    i=0
    for vec in dvecs:
        if vec == dvecs[omit]:
            continue
        for x in xrange(len(vec)):
            training[i] = vec[x]
            i = i + 1
    
    i=0
    for vec in cvecs:
        if vec == cvecs[omit]:
            continue 
        for x in xrange(len(vec)):
            classes[i] = vec[x]
            i = i + 1
            
    return training, classes


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def buildTrainingVector(data):
    dim = len(data["data"])
    training = DataVector(len(data["data"][0]), dim)
    
    for i in xrange(len(data["data"][0])):
        for d in xrange(dim):
            training[i*dim + d] = data["data"][d][i]
    
    return training

#-------------------------------------------------------------------------------
def buildYVector(data):
    y = DataVector(len(data["classes"]))
    for i in xrange(len(data["classes"])):
        y[i] = data["classes"][i]
        
    return y


#===============================================================================
# data_providers
#===============================================================================

def constructNormal(data):
    if len(data) != 1:
        raise Exception("Only one data file supported.")

    return data[0], len(data[0]["data"])

def splitNormal(data):
    yield (buildTrainingVector(data), buildYVector(data)), None
    return

## Random Fold
def constructFold(data):
    if len(data) != 1:
        raise Exception("Only one data file supported.")
    
    return split_n_folds(data[0], options.f_level, options.seed), len(data[0]["data"])

## Stratified Fold
def constructFoldr(data):
    if len(data) != 1:
        raise Exception("Only one data file supported.")
    
    return split_n_folds_stratified(data[0], options.f_level, options.seed), len(data[0]["data"])


## Sequential Fold
def constructFolds(data):
    if len(data) != 1:
        raise Exception("Only one data file supported.")
    
    return split_n_folds_sequential(data[0], options.f_level), len(data[0]["data"])


def splitFold(data):
    dvec, cvec = data
    for i in xrange(options.f_level):
        training = assembleTrainingVector(dvec,cvec,i)
        testing = (dvec[i], cvec[i])
        yield (training, testing)
    return

## list of data providers
# construct should bring the data in a suitable form.
# split should be a generator returning a tuple (training, testing) for the current fold.
# training/testing should be of the form (data, classes) or None if not present
data_providers = {
        "normal" : {"construct" : constructNormal, "split" : splitNormal},
        "fold" : {"construct" : constructFold, "split" : splitFold},
        "folds" : {"construct" : constructFolds, "split" : splitFold},                  
        "foldr" : {"construct" : constructFoldr, "split" : splitFold},                  
        }

#===============================================================================
# refine_providers
#===============================================================================

class SurplusRefineProvider:
    def __init__(self, status):
        self.status = status
        
    def reset(self):
        """Called before each learning step"""
        self.alpha = DataVector(self.status.grid.getStorage().size())
        self.alpha.setAll(0.0)
        self.i = 0
    
    def add(self, alpha, training):
        """Called once per fold"""
        self.alpha.add(alpha)
        self.i += 1
    
    def refine(self):
        """Called after learning all folds"""
        self.alpha.mult(1.0/self.i)
        
        refine = self.status.grid.createGridGenerator()
        
        functor = SurplusRefinementFunctor(self.alpha)
        
        refine.refine(functor)
        
        print "GridPoints: ", self.status.grid.getStorage().size()

##
# TODO: rework
class ErrorRefineProvider:
    def __init__(self, status):
        self.status = status
        
    def reset(self):
        """Called before each learning step"""
        self.alpha = DataVector(self.status.grid.getStorage().size())
        self.alpha.setAll(0.0)
        self.i = 0
    
    def add(self, alpha, training):
        """Called once per fold"""
        
        temp = DataVector(training[1].getSize())
        B = self.status.grid.createOperationB()
        
        B.multTranspose(alpha, training[0], temp)
        
        temp.sub(training[1])
        temp.sqr()        
        
        error_temp = DataVector(alpha.getSize())
        B.mult(temp, training[0], error_temp)
        
        self.alpha.add(error_temp)
        self.i += 1
    
    def refine(self):
        """Called after learning all folds"""
        self.alpha.mult(1.0/self.i)
        
        refine = self.status.grid.createGridGenerator()
        
        functor = SurplusRefinementFunctor(self.alpha)
        
        refine.refine(functor)
        
        print "GridPoints: ", self.status.grid.getStorage().size()

##list of refine providers
# see SuprlusRefineProvider for an example
refine_providers = {
        "surplus" : SurplusRefineProvider,
        "error" : ErrorRefineProvider,
        }

#===============================================================================
# eval_providers
#===============================================================================

def testVector(status, alpha, data, classes):
    eval = status.grid.createOperationEval()
    return eval.test(alpha, data, classes) / float(data.getSize())

class ClassesEvalProvider:
    def __init__(self, status):
        self.status = status
        
        self.training_overall = []
        self.testing_overall = []

    def reset(self):
        """Called before each learning step"""
        self.training_results = []
        self.testing_results = []

    def training(self, alpha, data):
        """Called to evaluate training data once per fold"""
        if data:
            self.training_results.append(testVector(self.status, alpha, data[0], data[1]))

    def testing(self, alpha, data):
        """Called to evaluate testing data once per fold"""
        if data:
            self.testing_results.append(testVector(self.status, alpha, data[0], data[1]))
                                     
    def evaluate(self):
        """Called after each learning step"""
        i = float(len(self.training_results))
        self.training_overall.append(sum(self.training_results)/i)
        self.testing_overall.append(sum(self.testing_results)/i)
        
        print self.training_overall
        print self.testing_overall

class RegressionEvalProvider:
    def __init__(self, status):
        self.status = status
        
        self.training_overall = []
        self.testing_overall = []

    def reset(self):
        """Called before each learning step"""
        self.training_results = []
        self.testing_results = []

    def training(self, alpha, data):
        """Called to evaluate training data once per fold"""
        if data:
            self.training_results.append(self.regressionTest(alpha, data[0], data[1]))

    def testing(self, alpha, data):
        """Called to evaluate testing data once per fold"""
        if data:
            self.testing_results.append(self.regressionTest(alpha, data[0], data[1]))
             
    def regressionTest(self, alpha, data, classes):
        temp = DataVector(classes.getSize())
        B = self.status.grid.createOperationB()
        B.multTranspose(alpha, data, temp)
        
        temp.sub(classes)
        temp.sqr()
        return temp.sum()/temp.getSize()
        
                                     
    def evaluate(self):
        """Called after each learning step"""
        i = float(len(self.training_results))
        self.training_overall.append(sum(self.training_results)/i)
        self.testing_overall.append(sum(self.testing_results)/i)
        
        print self.training_overall
        print self.testing_overall

## List of available refine providers
# See ClassesEvalProvider for an example        
eval_providers = {
        "classes" : ClassesEvalProvider,
        "regression" : RegressionEvalProvider,
        }

#===============================================================================
# General algorithms
#===============================================================================


class Status(object):
    """Grid status type. Contains everything needed for datamining"""
    def __init__(self, grid, mode, basetype, zeh):
        self.grid = grid
        self.mode = mode
        self.basetype = basetype
        self.zeh = zeh
        
        self.refine = None
        self.eval = None

    def __getstate__(self):
        odict = self.__dict__.copy()
        del odict['grid']
        
        del odict['refine']
        del odict['eval']
        
        odict['grid'] = self.grid.serialize()
        
        return odict

    def __setstate__(self, dict):
        #restore grid
        
        self.__dict__.update(dict)
        self.grid = Grid.unserialize(dict['grid'])

def learningStep(status, data):
    if status.eval:
        status.eval.reset()
    
    if status.refine:
        status.refine.reset()
        
    for training, testing in data_providers[status.mode]["split"](data):
        alpha = DataVector(status.grid.getStorage().size())
        alpha.setAll(0.0)

        m = Matrix(status.grid, training[0], options.l, status.zeh, status.basetype)
        b = m.generateb(training[1])

        res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, False, options.verbose)
        print res
        
        if status.eval:
            status.eval.training(alpha, training)
            status.eval.testing(alpha, testing)
        
        if status.refine:
            status.refine.add(alpha, training)
            
    if options.checkpoint:
        import pickle
        
        filename = "%s_%08d" % (options.checkpoint, status.grid.getStorage().size())
        fout = open(filename + ".grid", "wb")
        pickle.dump(status, fout, -1)
        fout.close()
        
        #writeAlphaARFF(filename + ".alpha.arff", alpha)


def learningAlgorithm(mode):
    data = loadData()
    data, dim = data_providers[mode]["construct"](data)

    status = construction(dim, mode)
    
    status.refine = refine_providers[options.refine](status)
    status.eval = eval_providers[options.eval](status)
    
    learningStep(status, data)
    if status.eval:
        status.eval.evaluate()
    
    for i in xrange(options.adaptive):
        if status.refine:
            status.refine.refine()
        learningStep(status, data)
        if status.eval:
            status.eval.evaluate()
    

def evalAlgorithm(mode):
    
    raise Exception("eval unsupported.")
    
    data = loadData()
    data, dim = data_providers[mode]["construct"](data)

    #TODO: construct grid from serialized file
    status = construction(dim, mode)

    alpha = DataVector(grid.size())
    #Fill alpha Vector
    
    for training, testing in data_providers[status.mode]["split"](data):
        status.eval.testing(alpha, testing)
        
    status.eval.evaluate()
    


#===============================================================================
# Main
#===============================================================================

#-------------------------------------------------------------------------------
def doUnsupported(mode):
    raise Exception(mode + " mode unsupported")

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
    parser.add_option("-m", "--mode", action="store", type="string", default="apply", dest="mode", help="Specifies the action to do. Get help for the mode please type --mode help.")
    parser.add_option("-C", "--zeh", action="store", type="string", default="laplaceadaptive", dest="zeh", help="Specifies the action to do.")
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
        print("Wrong C-mode! Please refer to --C help for further information.")
        sys.exit(1)

    # Execute the mode
    exec_mode(options.mode.lower())



if __name__=="__main__":
    _main()