############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2007-2009 Dirk Plueger (Dirk.Pflueger@in.tum.de)            #
# Copyright (C) 2007 Joerg Blank (blankj@in.tum.de)                         #
# Copyright (C) 2007 Richard Roettger (roettger@in.tum.de)                  #
# Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU Lesser General Public License as published  #
# by the Free Software Foundation; either version 3 of the License, or      #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU Lesser General Public License  #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################

## @package classifier
# @ingroup bin
# @brief Main script to do classification, regression, ...
# @version $CURR$
# @todo Join the different modes so that not 100 different versions (checkpointing, ...) exist!


from optparse import OptionParser
import sys
from tools import *
from pysgpp import *
from painlesscg import cg,sd,cg_new
from math import sqrt
from math import ceil
import random

from array import array

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
#
# @todo remove hack for level new when porting to classifier.new.py
#
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
                # hack for level 0
                if attrib2 == "level":
                    option = getattr(options, attrib2, None)
                    if option >= 0:
                        b = True
                else:
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
        print ("More on the usage of %s with --help" % (sys.argv[0]))
        sys.exit(1)

    # execute action
    modes[mode]['action']()
    


#-------------------------------------------------------------------------------
## Opens and read the data of an ARFF file
# Opens a file given by a filename
# @param filename filename of the file
# @return the data stored in the file as a set of arrays
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
# @todo Integrate into all modes
def constructGrid(dim):
    if options.grid == None:
        if options.uscaledboundary == True or options.completeboundary == True:
            if options.uscaledboundary == True:
                if options.verbose:
                    print "LinearBoundaryUScaledGrid, l=%s" % (options.level)
                grid = Grid.createLinearBoundaryUScaledGrid(dim)            
            if options.completeboundary == True:
                if options.verbose:
                    print "LinearBoundaryGrid, l=%s" % (options.level)
                grid = Grid.createLinearBoundaryGrid(dim)            
        else:
            if options.border:
                if options.polynom > 1:
                    if options.verbose:
                        print "ModPolyGrid, p=%d, l=%d" %(options.polynom, options.level)
                    grid = Grid.createModPolyGrid(dim, options.polynom)
                else:
                    if options.verbose:
                        print "ModLinearGrid, l=%s" % (options.level)
                    grid = Grid.createModLinearGrid(dim)
            else: #no border points
                if options.polynom > 1:
                    if options.verbose:
                        print "PolyGrid, p=%d, l=%d" %(options.polynom, options.level)
                    grid = Grid.createPolyGrid(dim, options.polynom)
                else:
                    if options.verbose:
                        print "LinearGrid, l=%s" % (options.level)
                    grid = Grid.createLinearGrid(dim)
                    
        generator = grid.createGridGenerator()
        generator.regular(options.level)
    else: #read grid from file
        if options.verbose:
            print "reading grid from %s" % (options.grid)
        grid = readGrid(options.grid)
    
    return grid


#-------------------------------------------------------------------------------
## Classify a dataset with an existing grid (compute accuracy).
# Evaluates a sparse grid function, given by a regular sparse grid (dim, level)
# and an alpha-Vector (ARFF-file) at a set of data points (ARFF-file).
# Outputs classifications to ARFF-file.
def doApply():
    # read data
    data = openFile(options.data[0])
    dim = len(data["data"])
    numData = len(data["data"][0])
    
    # read alpha vector
    alpha = buildTrainingVector(openFile(options.alpha))
    
    # construct corresponding grid
    grid = constructGrid(dim)
    if(alpha.getSize() != grid.getStorage().size()):
        print "Error: Inconsistent Grid and Alpha-Vector"
        print "alpha size %d, grid size %d" % (alpha.getSize(),grid.getStorage().size())
        sys.exit(1)
    
    # copy data to DataVector
    (x,y) = createDataVectorFromDataset(data)

    # evaluate
    q = DataVector(1,dim)
    classes = []
    # if test data contains classes, additionally compute accuracy
    if data.has_key("classes"):
        compute_accuracy = True
        acc = 0
    else:
        compute_accuracy = False
    # traverse Data
    for i in xrange(numData):
        x.getRow(i,q)
        val = grid.createOperationEval().eval(alpha,q)
        if compute_accuracy:
            if (val >= 0 and data["classes"][i] >= 0) or (val < 0 and data["classes"][i] < 0):
                acc += 1
        if val >= 0:
            classes.append(1.0)
        else:
            classes.append(-1.0)
    # output accuracy:
    acc = acc / float(numData)
    print "Accuracy on test data: %9.5f%%" % (100*acc)
    
    if options.regression: 
        m = Matrix(grid, buildTrainingVector(data), options.regparam, options.zeh)
        mse = evaluateError(buildYVector(data), alpha, m)[0]
            
    # get filename for output file
    if(options.outfile != None):
        data["filename"] = options.outfile
    else:
        data["filename"] = data["filename"] + ".out.arff"
    
    # write data with function values to ARFF-file
    data["classes"] = classes
    writeDataARFF([data])


#-------------------------------------------------------------------------------
## Evaluate a sparse grid function given by a grid and an 
# alpha-Vector (ARFF-file) at a set of data points (ARFF-file).
# If the data set contains class values, output some statistics 
# (classification accuracy, or - if regularization - mse and other measures).
# Outputs classification or (if regression) function values to ARFF-file.
def doEval():
    # read data
    data = openFile(options.data[0])
    dim = len(data["data"])
    numData = len(data["data"][0])

    # read alpha vector
    alpha = buildTrainingVector(openFile(options.alpha))

    # construct corresponding grid
    grid = constructGrid(dim)
    if(alpha.getSize() != grid.getStorage().size()):
        print "Error: Inconsistent Grid and Alpha-Vector"
        print "alpha size %d, grid size %d" % (alpha.getSize(),grid.getStorage().size())
        sys.exit(1)

    # copy data to DataVector
    (x,y) = createDataVectorFromDataset(data)

    # evaluate
    q = DataVector(1,dim)
    classes = []
    # if test data contains function values, additionally compute L2-norm of error
    if data.has_key("classes"):
        compute_accuracy = True
        err = 0
    else:
        compute_accuracy = False

    # classification:
    if not options.regression:
        # traverse Data
        for i in xrange(numData):
            x.getRow(i,q)
            val = grid.createOperationEval().eval(alpha,q)
            if compute_accuracy:
                if (val >= 0 and data["classes"][i] >= 0) or (val < 0 and data["classes"][i] < 0):
                    acc += 1
            if val >= 0:
                classes.append(1.0)
            else:
                classes.append(-1.0)
        if compute_accuracy:
            # output accuracy:
            acc = acc / float(numData)
            print "Accuracy on test data: %9.5f%%" % (100*acc)
    # regression:
    else: 
        # traverse Data
        for i in xrange(numData):
            x.getRow(i,q)
            val = grid.createOperationEval().eval(alpha,q)
            classes.append(val)
        if compute_accuracy:
            # output accuracy:
            m = Matrix(grid, x, options.regparam, options.zeh)
            evaluateError(y, alpha, m)[0]

    # get filename for output file
    if(options.outfile != None):
        data["filename"] = options.outfile
    else:
        data["filename"] = data["filename"] + ".out.arff"
    # write data with classes or function values to ARFF-file
    data["classes"] = classes
    writeDataARFF([data])

    
#-------------------------------------------------------------------------------
## Evaluate a sparse grid function given by grid and alphas at data points
# which are provided at stdin.
# Evaluates a sparse grid function, given by a regular sparse grid (dim and level
# or file) and an alpha-Vector (ARFF-file) at a set of data points (stdin).
# Data points have to be provided linewise via stdin (doubles).
# Outputs function values linewise to stdout.
def doEvalStdin():
    # read alpha vector
    alpha = buildTrainingVector(openFile(options.alpha))

    # construct corresponding grid
    if not (options.grid or (options.dim and options.level)):
        print "Error: No grid provided!"
        sys.exit(1)
    if options.dim:
        dim = options.dim
    else:
        dim = 0
    grid = constructGrid(dim)
    dim = grid.getDim()
    if(alpha.getSize() != grid.getStorage().size()):
        print "Error: Inconsistent Grid and Alpha-Vector"
        sys.exit(1)

    # read in normalization information if available
    if options.normfile:
        (border, minvals, maxvals, deltavals) = readNormfile(options.normfile)
        
    # evaluate
    q = DataVector(dim)
    # read data
    sys.stdout.write("---- enter data ---- stop with blank line ----\n")
    sys.stdout.flush()
    s = sys.stdin.readline().strip()
    while (s):
        dat = s.split(None)
        if len(dat) < dim:
            dat = s.split(",")
        for i in range(dim):
            if options.normfile:
                if deltavals[i] == 0:
                    q[i] = 0.5
                else:
                    q[i] = (float(dat[i])-minvals[i])/deltavals[i] + border
            else:
                q[i] = float(dat[i])
        val = grid.createOparationEval().eval(q, alpha)
        sys.stdout.write("%g\n" % val)
        sys.stdout.flush()
        s = sys.stdin.readline().strip()
    sys.stdout.write("\n")
    sys.stdout.flush()

    
#-------------------------------------------------------------------------------
## Learn a dataset
def doNormal():
    # read data
    data = openFile(options.data[0])
    dim = len(data["data"])
    numData = len(data["data"][0])
    
    if options.verbose:
        print "Dimension is:", dim
        print "Number of datasets is:", numData
        
    training = buildTrainingVector(data)
    y = buildYVector(data)

    grid = constructGrid(dim)

    alpha = run(grid, training, y)

#    for n in xrange(grid.getStorage().size()):
#        print grid.getStorage().get(n).getCoordinates()

    if options.outfile:
        writeAlphaARFF(options.outfile, alpha)
    
    if(options.gnuplot != None):
        if(dim != 2):
            print("Wrong dimension for gnuplot-Output!")
        else:
            writeGnuplot(options.gnuplot, grid, alpha, options.res)
        

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def run(grid, training, classes):
    alpha = None
    errors = None

    for adaptStep in xrange(options.adaptive + 1):
        print "Adaptive Step:", (options.adapt_start + adaptStep)
        alpha = DataVector(grid.getStorage().size())
        alpha.setAll(0.0)

        m = Matrix(grid, training, options.regparam, options.zeh)
        b = m.generateb(classes)
            
        res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, False, options.verbose, max_threshold=options.max_r)
        print "Conjugate Gradient output:"
        print res

        if options.regression:
            errors = evaluateError(classes, alpha, m)[1]

        if options.checkpoint != None:
            writeCheckpoint(options.checkpoint, grid, alpha, options.adapt_start + adaptStep)
        
        
        if(adaptStep  < options.adaptive):
            #if(options.verbose):
            print("refining grid")
            if options.regression:
                grid.createGridGenerator().refine(SurplusRefinementFunctor(errors, options.adapt_points))
            else:
                grid.createGridGenerator().refine(SurplusRefinementFunctor(alpha, options.adapt_points))
    return alpha

##
#Subroutine evaluation of error
#@todo remove printing messages from the subroutine and place it into the suited methods
#
def evaluateError(classes, alpha, m):   
    error = DataVector(len(classes))
    m.B.multTranspose(alpha, m.x, error)
    error.sub(classes) # error vector
    error.sqr() # entries squared
    # output some statistics
    err_min = error.min(1)
    err_max = error.max(1)
    print "(Min,Max) error: (%f,%f)" % (sqrt(err_min), sqrt(err_max))
    # output accuracy
    err = error.sum()
    print "L2-norm of error on data: %f" % (sqrt(err))
    mse = err / error.getSize()
    print "MSE: ", mse

    # calculate error per basis function
    errors = DataVector(alpha.getSize())
    m.B.mult(error, m.x, errors)
    
    return (mse, errors)


#-------------------------------------------------------------------------------
## Learn a dataset with a test dataset.
def doTest():
    data = openFile(options.data[0])
    test = openFile(options.test)

    training = buildTrainingVector(data)
    y = buildYVector(data)

    test_data = buildTrainingVector(test)
    test_classes = buildYVector(test)

    dim = len(data["data"])
    grid = constructGrid(dim)
    
    te_refine = []
    tr_refine = []
    num_refine = []
    
    adaptStep = 0
    while True: #loop exit condition on the end of the loop
        print "Adaptive Step:", (options.adapt_start + adaptStep)
        m = Matrix(grid, training, options.regparam, options.zeh)
        b = m.generateb(y)
        
        alpha = DataVector(grid.getStorage().size())
        alpha.setAll(0.0)
        res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, False, options.verbose, max_threshold=options.max_r)
        print "Conjugate Gradient output:"
        print res

        tr = testVectorFast(grid, alpha, training, y)
        te = testVectorFast(grid, alpha, test_data, test_classes)
        num_refine.append(grid.getStorage().size())
        
        if options.regression:
            print "Training Data:"
            tr, errors = evaluateError(y, alpha, m)
            tr_refine.append(tr)
            print "Test Data:"
            m_test = Matrix(grid, test_data, options.regparam, options.zeh)
            mse = evaluateError(test_classes, alpha, m_test)[0]
            te_refine.append(mse)
        else:     
            tr_refine.append(tr)
            te_refine.append(te)
        
        if options.verbose and not options.regression:
            print "Correct classified on training data: ",tr
            print "Correct classified on testing data:  ",te

        if options.checkpoint != None: writeCheckpoint(options.checkpoint, grid, alpha, (options.adapt_start + adaptStep))
        if options.stats != None: writeStats(options.stats, formTxt([te_refine[-1]], [tr_refine[-1]], [num_refine[-1]], False))
            
        #increment adaptive step
        adaptStep += 1
        if(options.adaptive >= 0 and adaptStep <= options.adaptive) \
            or (options.epochs_limit > 0 and options.epochs_limit > getEpochsErrorIncreasing(te_refine)) \
            or (options.regression and options.mse_limit > 0 and options.mse_limit > te_refine[-1]) \
            or (options.grid_limit > 0 and options.grid_limit > grid.getStorage().size()) \
            :
            print("refining grid")
            
            numOfPoints = 0
            if options.adapt_rate: 
                numOfPoints = int(ceil( options.adapt_rate * grid.createGridGenerator().getNumberOfRefinablePoints()))
            else: numOfPoints = options.adapt_points
            
            if options.regression:
                grid.createGridGenerator().refine(SurplusRefinementFunctor(errors, numOfPoints))
            else:
                grid.createGridGenerator().refine(SurplusRefinementFunctor(alpha, numOfPoints))
           
            if(options.verbose): print("Number of points: %d" %(grid.getStorage().size(),))
            
        #Break condition for while-loop: 
        #number of adaptive steps is achieved or MSE of test data increase last 20 epochs or last MSE less then given boundary 
        else:
            break
        
    #--end of while loop

    if options.stats != None:
        txt = formTxt(te_refine, tr_refine, num_refine)
        #writeStats(options.stats, txt)
        if options.verbose: print txt
    
    return (tr_refine, te_refine, num_refine)

##
# returns the number of epochs the error is increasing
# @param list: List with MSE's from different refinement iterations
#
def getEpochsErrorIncreasing(list):
    length = len(list)
    if length == 0:
        return 0
    else:
        for i in xrange(1,length):
            if list[-i] < list[-i - 1]:
                return i
    return length

def formHeader():
    return  "#\t%f, %-10g, %f" % (options.level, options.regparam, options.adaptive)

###
## returns txt variable for stats and checkpoint
##
def formTxt(te_refine, tr_refine, num_refine, withHeader = True):
    txt = ""
    #if withHeader:
    #    txt = formHeader()
    for i in xrange(len(tr_refine)):
        txt = txt + "%f %.10f %f %f %f %.10f %.10f" % (options.level, options.regparam, options.adaptive, i, num_refine[i], tr_refine[i], te_refine[i])
    
    return txt + "\n"

##
# returns txt variable for stats
#
def formTxtVal(te_refine, val_refine, tr_refine, num_points, withHeader = True):
    txt = "%d, %-10g, %f" % (options.level, options.regparam, options.adaptive)
    for i in xrange(len(tr_refine)):
        txt = txt + ", %f, %.10f, %.10f, %.10f" % (num_points[i], tr_refine[i], val_refine[i], te_refine[i])
    return txt + "\n"

#-------------------------------------------------------------------------------
## Learn a dataset with a random n-fold.
def doFold():
    data = openFile(options.data[0])
    (dvec, cvec) = split_n_folds(data, options.f_level, options.seed)
        
    if options.regression:
        performFoldRegression(dvec, cvec)
    else:
        performFold(dvec, cvec)


#-------------------------------------------------------------------------------
## Learn a dataset with a sequential n-fold.
def doFolds():
    data = openFile(options.data[0])
    (dvec, cvec) = split_n_folds_sequential(data, options.f_level)

    if options.regression:
        performFoldRegression(dvec, cvec)
    else:
        performFold(dvec, cvec)
    return

#-------------------------------------------------------------------------------
## Learn a dataset with a stratified n-fold.
def doFoldr():
    data = openFile(options.data[0])
    (dvec, cvec) = split_n_folds_stratified(data, options.f_level, options.seed)
    
    if options.regression:
        performFoldRegression(dvec, cvec)
    else:
        performFold(dvec, cvec)

#-------------------------------------------------------------------------------
## Learn a dataset with a stratified n-fold.
def doFoldStratified():
    data = openFile(options.data[0])
    (dvec, cvec) = split_n_folds_stratified(data, options.f_level, options.seed)

    if options.regression:
        raise Exception("Not implemented!")
    else:
        # run all folds?
        if not options.onlyfoldnum:
            for i in range(options.f_level):
                performFoldNew(dvec, cvec, i)
        # run only single one of the folds
        else:
            performFoldNew(dvec, cvec, options.onlyfoldnum)


#-------------------------------------------------------------------------------
## Learn a dataset with a n-fold from a set of files.
def doFoldf():

    if len(options.data)==1:
        print("Error: Not enough to do. At least two data-Files necessary!")
        sys.exit(1)
        
    dvec = []
    cvec = []
    
    for file in options.data:
        data = openFile(file)
        dvec.append(buildTrainingVector(data))
        cvec.append(buildYVector(data))
        
    options.f_level = len(dvec)
    
    if options.regression:
        performFoldRegression(dvec, cvec)
    else:
        performFold(dvec, cvec)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def performFold(dvec,cvec):

    grid = constructGrid(dvec[0].getDim())
        
    num_points = []
    tr_refine = []
    te_refine = []
    
    for adaptStep in xrange(options.adaptive + 1):
        trainingCorrect = []
        testingCorrect =[]

        refinealpha = DataVector(grid.getStorage().size())
        refinealpha.setAll(0.0)

        alpha = DataVector(refinealpha)
        
        for foldSetNumber in xrange(options.f_level):
#            alpha.setAll(0.0)
            training,classes = assembleTrainingVector(dvec,cvec,foldSetNumber)
            
            m = Matrix(grid, training, options.regparam, options.zeh)
            b = m.generateb(classes)

            res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, options.reuse, options.verbose, max_threshold=options.max_r)
            print res

            tr = testVectorFast(grid, alpha, training, classes)
            te = testVectorFast(grid, alpha, dvec[foldSetNumber], cvec[foldSetNumber])

            trainingCorrect.append(tr)
            testingCorrect.append(te)

# Verstehe ich nicht. Habs deshalb auskommentiert (Dirk)
#             if(adaptStep +1 == options.adaptive):
#                 #Letzte verfeinerung, wir sind fertig
#                 pass
#             else:
#                 refinealpha.add(alpha)
            refinealpha.add(alpha)

        if options.verbose:
            print(trainingCorrect)   
            print(testingCorrect)

      
        tr = sum(trainingCorrect)/options.f_level
        te = sum(testingCorrect)/options.f_level
        
        if options.verbose:
            print "training: ",tr
            print "testing:  ",te

        num_points.append(grid.getStorage().size())
        tr_refine.append(tr)
        te_refine.append(te)
        
        refinealpha.mult(1.0/options.f_level)
        
        if options.checkpoint != None: writeCheckpoint(options.checkpoint, grid, refinealpha)
        
        if options.stats != None: writeStats(options.stats, formTxt(te_refine, tr_refine, num_points)) 
        
        if options.verbose:
            print "alpha"
            print refinealpha
            print "grid"
            print grid
        if(adaptStep < options.adaptive):
            print "refine"
            grid.createGridGenerator().refine(SurplusRefinementFunctor(refinealpha))

    if options.stats != None:
        txt = formTxt(te_refine, tr_refine, num_points)
        writeStats(options.stats, txt)
        if options.verbose: print txt

    return

#-------------------------------------------------------------------------------
## Perform n-fold cross validation for fold ifold.
# Splits the training data set into a training and a validation set
# @param dvec List of DataVectors, containing the test data for each fold
# @param cvec List of DataVectors, containing the classes of the test data for each fold
# @param ifold Number of fold to perform
def performFoldNew(dvec,cvec,ifold):

    print "Starting fold %d" % (ifold)

    # init
    grid = constructGrid(dvec[0].getDim())
    num_points = []
    tr_refine = []
    val_refine = []
    te_refine = []
    
    # loop until number of refinements reached
    for adaptStep in xrange(options.adaptive + 1):
        if options.verbose: print "Step %d" % (adaptStep)
        trainingCorrect = []
        testingCorrect =[]

        # construct/split DataVectors
        training,classes = assembleTrainingVector(dvec, cvec, ifold)
        data_tr,data_val,class_tr,class_val = split_DataVectors_by_proportion_stratified(training, classes, 0.66)

        # construct and solve CG
        alpha = DataVector(grid.getStorage().size())
        alpha.setAll(0.0)
        m = Matrix(grid, data_tr, options.regparam, options.zeh)
        b = m.generateb(class_tr)
        res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, options.reuse, options.verbose, max_threshold=options.max_r)
        if options.verbose: print res

        # training and validation accuracy
        tr = testVectorFast(grid, alpha, data_tr, class_tr)
        val = testVectorFast(grid, alpha, data_val, class_val)

        # compute accuracy on test set
        # Therefore construct and solve CG again for whole training data and evaluate on test data
        m = Matrix(grid, training, options.regparam, options.zeh)
        b = m.generateb(classes)
        res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, True, options.verbose, max_threshold=options.max_r)
        te = testVectorFast(grid, alpha, dvec[ifold], cvec[ifold])

        num_points.append(grid.getStorage().size())
        tr_refine.append(tr)
        val_refine.append(val)
        te_refine.append(te)

        if options.verbose:
            print "num_points:", grid.getStorage().size()
            print "training:  ", tr
            print "validating:", val
            print "testing:   ", te

        if options.checkpoint != None: writeCheckpoint(options.checkpoint, grid, alpha, options.adapt_start + adaptStep, fold=ifold)
        
        if(adaptStep < options.adaptive):
            if options.verbose: print "refining"
            grid.createGridGenerator().refine(SurplusRefinementFunctor(alpha))

    txt = formTxtVal(tr_refine, val_refine, te_refine, num_points)
    print txt
    if options.stats != None:
        writeStats(options.stats, txt)

    return

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def performFoldRegression(dvec,cvec):
    
    grid = constructGrid(dvec[0].getDim())
        
    num_points = []
    tr_refine = []
    te_refine = []
    tr_meanSqrError = []
    te_meanSqrError = []
        
    for adpatStep in xrange(options.adaptive + 1):
        trainingCorrect = []
        testingCorrect = []
        meanSqrErrorsTraining = []
        meanSqrErrorsTesting = []        

        refineerrors = DataVector(grid.getStorage().size())
        refineerrors.setAll(0.0)

        alpha = DataVector(grid.getStorage().size())
        
        for foldSetNumber in xrange(options.f_level):
#            alpha.setAll(0.0)
            training,classes = assembleTrainingVector(dvec,cvec,foldSetNumber)
            
            m = Matrix(grid, training, options.regparam, options.zeh)
            b = m.generateb(classes)

            res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, options.reuse, options.verbose, max_threshold=options.max_r)
            print res
            tr = testVector(grid,alpha,training,classes)
            te = testVector(grid,alpha,dvec[foldSetNumber],cvec[foldSetNumber])

            trainingCorrect.append(tr)
            testingCorrect.append(te)
            
            # calculate Mean Square Error for training set
            temp = DataVector(classes.getSize())
            m.B.multTranspose(alpha, m.x, temp)
            temp.sub(classes)
            temp.sqr()
            meanSqrErrorsTraining.append(temp.sum() / temp.getSize())

            # calculate error per base function
            errors = DataVector(alpha.getSize())
            m.B.mult(temp, m.x, errors)

            # calculate Mean Square Error for testing set
            temp = DataVector(cvec[foldSetNumber].getSize())
            m.B.multTranspose(alpha, dvec[foldSetNumber], temp)
            temp.sub(cvec[foldSetNumber])
            temp.sqr()
            meanSqrErrorsTesting.append(temp.sum() / temp.getSize())

            
            if(adpatStep +1 == options.adaptive):
                #Letzte verfeinerung, wir sind fertig
                pass
            else:
                refineerrors.add(errors)

        if options.verbose:
            print(trainingCorrect)   
            print(testingCorrect)
            print(meanSqrErrorsTraining)
            print(meanSqrErrorsTesting)
            
      
        tr = sum(trainingCorrect)/options.f_level
        te = sum(testingCorrect)/options.f_level
        trSqrError = sum(meanSqrErrorsTraining)/options.f_level
        teSqrError = sum(meanSqrErrorsTesting)/options.f_level
        
        if options.verbose:
            print "training: ",tr, trSqrError
            print "testing:  ",te, teSqrError

        num_points.append(grid.getStorage().size())
        tr_refine.append(tr)
        te_refine.append(te)
        tr_meanSqrError.append(trSqrError)
        te_meanSqrError.append(teSqrError)
        
        refineerrors.mult(1.0/options.f_level)
        
        if(adpatStep + 1 < options.adaptive):
            print "refine"
            grid.createGridGenerator().refine(SurplusRefinementFunctor(refineerrors))

    if options.stats != None:
            txt = formTxt(te_meanSqrError, tr_meanSqrError, num_points)
            writeStats(options.stats, txt )
            if options.verbose: print txt

    return

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
## Builds the training data vector
# 
# @param data a list of lists that contains the points a the training data set, coordinate-wise
# @return a instance of a DataVector that stores the training data
def buildTrainingVector(data):
    dim = len(data["data"])
    training = DataVector(len(data["data"][0]), dim)
    
    # i iterates over the data points, d over the dimension of one data point
    for i in xrange(len(data["data"][0])):
        for d in xrange(dim):
            training[i*dim + d] = data["data"][d][i]
    
    return training

#-------------------------------------------------------------------------------
## Computes the classification accuracy on some test data. 
# Tests on the classes {+1, -1}, cut-off at 0.
# @param grid the sparse grid
# @param alpha DataVector of surplusses
# @param test a DataVector containing a dataset of points to test on
# @param classes DataVector of correct class values
# @return classification accuracy
def testVector(grid,alpha,test,classes):
    p = DataVector(1,test.getDim())
    correct = 0
    for i in xrange(test.getSize()):
        test.getRow(i,p)
        val = grid.createOperationEval().eval(alpha,p)
        if (val < 0 and classes[i] < 0 ) or (val > 0 and classes[i] > 0 ):
            correct = correct + 1
            
    return float(correct)/test.getSize()

#-------------------------------------------------------------------------------
## Computes the classification accuracy on some test data. 
#
# Tests on the classes {+1, -1}, cut-off at 0. testVectorFast uses an OpenMP enabled c++ routine for testing
# @param grid the sparse grid
# @param alpha DataVector of surplusses
# @param test a DataVector containing a dataset of points to test on
# @param classes DataVector of correct class values
# @return classification accuracy
def testVectorFast(grid, alpha, test, classes):
    return grid.createOperationEval().test(alpha, test, classes)/float(test.getSize())

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def testVectorValues(grid,alpha,test,classes,evalValues):
    p = DataVector(1,test.getDim())
    correct = 0
    for i in xrange(test.getSize()):
        test.getRow(i,p)
        val = grid.EvaluatePoint(p, alpha)
        evalValues.append(val)
        if (val < 0 and classes[i] < 0 ) or (val > 0 and classes[i] > 0 ):
            correct = correct + 1
            
    return float(correct)/test.getSize()

#-------------------------------------------------------------------------------
## builds a vector that contains the class information
#
# @param data class information
# @return DataVector that contains the class information, length of the vector is identical to the length of the data vector
def buildYVector(data):
    y = DataVector(len(data["classes"]))
    for i in xrange(len(data["classes"])):
        y[i] = data["classes"][i]
        
    return y



#===============================================================================
# Main
#===============================================================================

# check so that file can also be imported in other files
if __name__=='__main__':

    # Initialize OptionParser, set Options
    parser = OptionParser()
    parser.add_option("-l", "--level", action="store", type="int", dest="level", help="Gridlevel")
    parser.add_option("-D", "--dim", action="callback", type="int",dest="dim", help="Griddimension", callback=callback_deprecated)
    parser.add_option("-a", "--adaptive", action="store", type="int", default=0, dest="adaptive", metavar="NUM", help="Using an adaptive Grid with NUM of refines")
    parser.add_option("--adapt_points", action="store", type="int", default=1, dest="adapt_points", metavar="NUM", help="Number of points in one refinement iteration")
    parser.add_option("--adapt_rate", action="store", type="int", dest="adapt_rate", metavar="NUM", help="Percentage of points from all refinable points in one refinement iteration")
    parser.add_option("--adapt_start", action="store", type="int", default=0, dest="adapt_start", metavar="NUM", help="The index of adapt step to begin with")
    parser.add_option("-m", "--mode", action="store", type="string", default="apply", dest="mode", help="Specifies the action to do. Get help for the mode please type --mode help.")
    parser.add_option("-C", "--zeh", action="store", type="string", default="laplace", dest="zeh", help="Specifies the action to do.")
    parser.add_option("-f", "--foldlevel", action="store", type="int",default=10, metavar="LEVEL", dest="f_level", help="If a fold mode is selected, this specifies the number of sets generated")
    parser.add_option("--onlyfoldnum", action="store", type="int", default=None, metavar="I", dest="onlyfoldnum", help="Run only fold I in n-fold cross-validation")
    parser.add_option("-L", "--lambda", action="store", type="float",default=0.000001, metavar="LAMBDA", dest="regparam", help="Lambda")
    parser.add_option("-i", "--imax", action="store", type="int",default=500, metavar="MAX", dest="imax", help="Max number of iterations")
    parser.add_option("-r", "--accuracy", action="store", type="float",default=0.0001, metavar="ACCURACY", dest="r", help="Specifies the accuracy of the CG-Iteration")
    parser.add_option("--max_accuracy", action="store", type="float", default=None, metavar="ACCURACY", dest="max_r", help="If the norm of the residuum falls below ACCURACY, stop the CG iterations")
    parser.add_option("-d", "--data", action="append", type="string", dest="data", help="Filename for the Datafile.")
    parser.add_option("-t", "--test", action="store", type="string", dest="test", help="File containing the testdata")
    parser.add_option("-A", "--alpha", action="store", type="string", dest="alpha", help="Filename for a file containing an alpha-Vector")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="Filename where the calculated alphas are stored")
    parser.add_option("-g", "--gnuplot", action="store", type="string", dest="gnuplot", help="In 2D case, the generated can be stored in a gnuplot readable format.")
    parser.add_option("-R", "--resolution", action="store", type="int",default=50, metavar="RESOLUTION", dest="res", help="Specifies the resolution of the gnuplotfile")
    parser.add_option("-s", "--stats", action="store", type="string", dest="stats", help="In this file the statistics from the test are stored")
    parser.add_option("-p", "--polynom", action="store", type="int", default=0, dest="polynom", help="Sets the maximum degree for high order basis functions. Set to 2 or larger to activate. Works only with 'identity' and 'fold'-modes.")
    parser.add_option("-b", "--border", action="store_true", default=False, dest="border", help="Enables special border base functions")
    parser.add_option("--uscaled-boundary", action="store_true", default=False, dest="uscaledboundary", help="Enables boundary functions that have a point on the boundary for every inner point")
    parser.add_option("--complete-boundary", action="store_true", default=False, dest="completeboundary", help="Enables boundary functions that have more points on the boundary than inner points")
    parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose", help="Provides extra output")
    parser.add_option("--normfile", action="store", type="string", dest="normfile", metavar="FILE", help="For all modes that read data via stdin. Normalizes data according to boundaries in FILE")
    parser.add_option("--reuse", action="store_true", default=False, dest="reuse", help="Reuse alpha-values for CG")
    parser.add_option("--seed", action="store", type="float", dest="seed", help="Random seed used for initializing")
    parser.add_option("--regression", action="store_true", default=False, dest="regression", help="Use regression approach.")
    parser.add_option("--checkpoint", action="store", type="string", dest="checkpoint", help="Filename for checkpointing. For fold? and test. No file extension.")
    parser.add_option("--grid", action="store", type="string", dest="grid", help="Filename for Grid-resume. For fold? and test. Full filename.")
#TODO: maybe a better name for parameter
    parser.add_option("--epochs_limit", action="store", type="int", default="0", dest="epochs_limit", help="Number of refinement iterations (epochs), MSE of test data have to increase, before refinement will stop.")
    parser.add_option("--mse_limit", action="store", type="float", default="0.0", dest="mse_limit", help="If MSE of test data fall below this limit, refinement will stop.")
    parser.add_option("--grid_limit", action="store", type="int", default="0", dest="grid_limit", help="If the number of points on grid exceed grid_limit, refinement will stop.")

    # parse options
    (options,args)=parser.parse_args()

    # check some options
    zeh = options.zeh.lower()

    #this incrementation is unobvious, the default value of options.adaptive 
    #options.adaptive = options.adaptive + 1

    # check C-mode
    if zeh == "help":
        print "The following C-modes are available:"
        for m in zeh_modes.keys():
            print "%15s: %s" % (m, zeh_modes[m])
        sys.exit(0)
    elif zeh not in zeh_modes.keys():
        print("Wrong C-mode! Please refer to --C help for further information.")
        sys.exit(1)

    if options.polynom > 1 and zeh != "identity":
        print("Wrong C-mode selected for high-order grids.")
        sys.exit(1)

    if options.polynom > 1 and options.border:
        print("Special border bases do not work for High-Order Grids")
        sys.exit(1)

    # check further parameters:
    if options.onlyfoldnum and (not options.onlyfoldnum in range(options.f_level)):
        raise Exception("--onlyfoldnum: Not in range 0,...,--foldlevel")



    # specifiy the modes:
    # modes is an array containing all modes, the options needed by the mode and the action
    # that is to be executed
    modes = {
        'apply'    : {'help': "classify a dataset with an existing grid (compute accuracy)",
                      'required_options': ['data', 'alpha', ['level', 'grid']],
                      'action': doApply},
        'eval'     : {'help': "evaluate a sparse grid function given by grid and alphas at data points",
                      'required_options': ['data', 'alpha', ['level', 'grid']],
                      'action': doEval},
        'evalstdin': {'help': "evaluate a sparse grid function given by grid and alphas at data points",
                      'required_options': ['alpha', ['level', 'grid']],
                      'action': doEvalStdin},
        'normal'   : {'help': "learn a dataset",
                      'required_options': ['data', ['level', 'grid']],
                      'action': doNormal},
        'test'     : {'help': "learn a dataset with a test dataset",
                      'required_options': ['data', 'test', ['level', 'grid']],
                      'action': doTest},
        'fold'     : {'help': "learn a dataset with a random n-fold",
                      'required_options': ['data', ['level', 'grid']],
                      'action': doFold},
        'folds'    : {'help': "learn a dataset with a sequential n-fold",
                      'required_options': ['data', ['level', 'grid']],
                      'action': doFolds},
        'foldstratified'    : {'help': "learn a dataset with a stratified n-fold",
                      'required_options': ['data', ['level', 'grid']],
                      'action': doFoldStratified},
        'foldr'    : {'help': "learn a dataset with a stratified n-fold",
                      'required_options': ['data', ['level', 'grid']],
                      'action': doFoldr},
        'foldf'    : {'help': "learn a dataset with a n-fold from a set of files",
                      'required_options': ['data', ['level', 'grid']],
                      'action': doFoldf}
        }

    # Execute the mode
    exec_mode(options.mode.lower())
