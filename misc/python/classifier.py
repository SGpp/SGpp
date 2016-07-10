#!/usr/bin/python

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

## @package classifier
# @ingroup bin
# @author Dirk Pflueger, Joerg Blank, Richard Roettger
# @brief Main script to do classification, regression, ...
# @version $CURR$
# @todo Join the different modes so that not 100 different versions (checkpointing, ...) exist!


from optparse import OptionParser
import sys,os

from pysgpp.extensions.datadriven.tools import *
from pysgpp import *
from pysgpp.extensions.datadriven import *
from painlesscg import cg,sd,cg_new
from math import sqrt
from math import ceil
import random

from array import array

try:
    import psyco
    psyco.full()
    print "#Using psyco"
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
## Opens and read the alphas of an ARFF (or plain whitespace-separated data) file.
# Opens a file given by a filename.
# @param filename filename of the file
# @return the data stored in the file as a set of arrays
def openAlphaFile(filename):
    return readAlpha(filename)


#-------------------------------------------------------------------------------
## Opens and read the data of an ARFF (or plain whitespace-separated data) file.
# Opens a file given by a filename.
# @param filename filename of the file
# @return the data stored in the file as a set of arrays
def openFile(filename):
    return readData(filename)

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
        # grid points on boundary
        if options.trapezoidboundary == True or options.completeboundary == True:
            if options.polynom > 1:
                print "Error. Not implemented yet."
                sys.exit(1)
            if options.trapezoidboundary == True:
                if options.verbose:
                    print "LinearBoundaryGrid, l=%s" % (options.level)
                grid = Grid.createLinearBoundaryGrid(dim)
            if options.completeboundary == True:
                if options.verbose:
                    print "LinearL0BoundaryGrid, l=%s" % (options.level)
                grid = Grid.createLinearBoundaryGrid(dim, 0)
        elif options.function_type == "modWavelet":
            if options.verbose:
                print "ModWaveletGrid, l=%s" % (options.level)
            grid = Grid.createModWaveletGrid(dim)
        else:
            # modified boundary functions?
            if options.border:
                if options.polynom > 1:
                    if options.verbose:
                        print "ModPolyGrid, p=%d, l=%d" %(options.polynom, options.level)
                    grid = Grid.createModPolyGrid(dim, options.polynom)
                else:
                    if options.verbose:
                        print "ModLinearGrid, l=%s" % (options.level)
                    grid = Grid.createModLinearGrid(dim)
            # grid points on boundary?
            elif options.boundary == 1:
                if options.polynom > 1:
                    print "Error. Not implemented yet."
                    sys.exit(1)
                else:
                    if options.verbose:
                        print "LinearBoundaryGrid, l=%s" % (options.level)
                    grid = Grid.createLinearBoundaryGrid(dim)
            # more grid points on boundary?
            elif options.boundary == 2:
                if options.polynom > 1:
                    print "Error. Not implemented yet."
                    sys.exit(1)
                else:
                    if options.verbose:
                        print "LinearL0BoundaryGrid, l=%s" % (options.level)
                    grid = Grid.createLinearBoundaryGrid(dim, 0)
            else: #no border points
                if options.polynom > 1:
                    if options.verbose:
                        print "PolyGrid, p=%d, l=%d" %(options.polynom, options.level)
                    grid = Grid.createPolyGrid(dim, options.polynom)
                else:
                    if options.verbose:
                        print "LinearGrid, l=%s" % (options.level)
                    grid = Grid.createLinearGrid(dim)

        generator = grid.getGenerator()
        generator.regular(options.level)
    else: #read grid from file
        if options.verbose:
            print "reading grid from %s" % (options.grid)
        grid = readGrid(options.grid)

    return grid

## Calculates the number of points, that should be refined
# @param options: options object
# @param grid: grid
# @return: number of points, that should be refined
def getNumOfPoints(options, grid):
    numOfPoints = 0
    if options.adapt_rate:
        numOfPoints = int(ceil( options.adapt_rate * grid.getGenerator().getNumberOfRefinablePoints()))
    else: numOfPoints = options.adapt_points
    return numOfPoints


#-------------------------------------------------------------------------------
## Classify a dataset with an existing grid (compute accuracy).
# Evaluates a sparse grid function, given by a regular sparse grid (dim, level)
# and an alpha-Vector (ARFF-file) at a set of data points (ARFF-file).
# Outputs classifications to ARFF-file.
def doApply():
    # read data
    data = openFile(options.data[0])
    dim = data["data"].getNcols()
    numData = data["data"].getNrows()

    # read alpha vector
#    alpha = buildTrainingVector(openAlphaFile(options.alpha))
    alpha = openAlphaFile(options.alpha)

    # construct corresponding grid
    grid = constructGrid(dim)
    if(len(alpha) != grid.getSize()):
        print "Error: Inconsistent Grid and Alpha-Vector"
        print "alpha size %d, grid size %d" % (len(alpha),grid.getSize())
        sys.exit(1)

    # copy data to DataVector
    #(x,y) = createDataVectorFromDataset(data)
    x = data['data']
    y = data['classes']

    # evaluate
    q = DataVector(dim)
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
        val = createOperationEval(grid).eval(alpha,q)
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
        m = Matrix(grid, buildTrainingVector(data), options.regparam, options.CMode, options.Hk)
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
    dim = data["data"].getNcols()
    numData = data["data"].getNrows()

    # read alpha vector
#    alpha = buildTrainingVector(openAlphaFile(options.alpha))
    alpha = openAlphaFile(options.alpha)

    # construct corresponding grid
    grid = constructGrid(dim)
    if(len(alpha) != grid.getSize()):
        print "Error: Inconsistent Grid and Alpha-Vector"
        print "alpha size %d, grid size %d" % (len(alpha),grid.getSize())
        sys.exit(1)

    # copy data to DataVector
#    (x,y) = createDataVectorFromDataset(data)
    x = data['data']
    y = data['classes']

    # evaluate
    q = DataVector(dim)
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
        acc = 0
        for i in xrange(numData):
            x.getRow(i,q)
            val = createOperationEval(grid).eval(alpha,q)
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
            val = createOperationEval(grid).eval(alpha,q)
            classes.append(val)
        if compute_accuracy:
            # output accuracy:
            m = Matrix(grid, x, options.regparam, options.CMode, options.Hk)
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
#    alpha = buildTrainingVector(openAlphaFile(options.alpha))
    alpha = openAlphaFile(options.alpha)

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
    if(len(alpha) != grid.getSize()):
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
    dim = data["data"].getNcols()
    numData = data["data"].getNrows()

    if options.verbose:
        print "Dimension is:", dim
        print "Size of datasets is:", numData
        print "Gridsize is:", grid.getSize()

    training = buildTrainingVector(data)
    y = buildYVector(data)

    grid = constructGrid(dim)

    alpha = run(grid, training, y)

#    for n in xrange(grid.getSize()):
#        print grid.getStorage().getPoint(n).getCoordinates()

    if options.outfile:
        writeAlphaARFF(options.outfile, alpha)
    if options.gridfile:
        writeGrid(options.gridfile, grid)

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
        alpha = DataVector(grid.getSize())
        alpha.setAll(0.0)

        m = Matrix(grid, training, options.regparam, options.CMode, options.Hk)
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
                grid.getGenerator().refine(SurplusRefinementFunctor(errors, getNumOfPoints(options, grid)))
            else:
                grid.getGenerator().refine(SurplusRefinementFunctor(alpha, getNumOfPoints(options, grid)))
    return alpha

##
#Subroutine evaluation of error
#@todo remove printing messages from the subroutine and place it into the suited methods
#
def evaluateError(classes, alpha, m):
    error = DataVector(len(classes))
    m.B.mult(alpha, error)
    error.sub(classes) # error vector
    error.sqr() # entries squared
    # output some statistics
    err_min = error.min()
    err_max = error.max()
    print "(Min,Max) error (abs.): (%f,%f)" % (sqrt(err_min), sqrt(err_max))
    # output accuracy
    mse = error.sum() / float(len(error))
    print "MSE: %g on %d data pts" % (mse, len(classes))
    print "RMSE / L2-norm of error on data: %g" % (sqrt(mse))

    # output functional
    N = len(alpha)
    temp = DataVector(N)
    m.C.mult(alpha, temp)
    Cnorm = alpha.dotProduct(temp)
    M = len(classes)
    temp2 = DataVector(M)
    m.B.mult(alpha, temp2)
    m.B.multTranspose(temp2, temp)
    BBTnorm = alpha.dotProduct(temp)
    print "functional: %g + %g * %g = %g" % (mse, m.l, Cnorm, mse+m.l*Cnorm)

    # calculate error per basis function
    errors = DataVector(len(alpha))
    m.B.multTranspose(error, errors)
    errors.componentwise_mult(alpha)
    return (mse, errors)


#-------------------------------------------------------------------------------
## Learn a dataset with a test dataset.
def doTest():
    if options.verbose:
        print "doTest here"
        print "data: %s\ntest: %s" % (options.data[0],options.test)
    data = openFile(options.data[0])
    test = openFile(options.test)

    training = buildTrainingVector(data)
    y = buildYVector(data)
    test_data = buildTrainingVector(test)
    test_classes = buildYVector(test)

    dim = data["data"].getNcols()
    grid = constructGrid(dim)
    if options.verbose:
        print "Dimension is:", dim
        print "Size of datasets is:", training.getNrows()
        print "Size of test datasets is:", test_data.getNrows()
        print "Gridsize is:", grid.getSize()

    te_refine = []
    tr_refine = []
    num_refine = []

    adaptStep = 0
    while True: #loop exit condition on the end of the loop
        print "Adaptive Step:", (options.adapt_start + adaptStep)
        m = Matrix(grid, training, options.regparam, options.CMode, options.Hk)
        b = m.generateb(y)

        alpha = DataVector(grid.getSize())
        alpha.setAll(0.0)
        res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, False, options.verbose, max_threshold=options.max_r)
        print "Conjugate Gradient output:"
        print res

        num_refine.append(grid.getSize())

        if options.regression:
            print "Training Data:"
            tr, errors = evaluateError(y, alpha, m)
            tr_refine.append(tr)
            print "Test Data:"
            m_test = Matrix(grid, test_data, options.regparam, options.CMode, options.Hk)
            mse = evaluateError(test_classes, alpha, m_test)[0]
            te_refine.append(mse)
        else:
            print "Training results:"
            tr = testVectorFastWithCharacteristicNumbers(grid, alpha, training, y)
            print "Test results:"
            te = testVectorFastWithCharacteristicNumbers(grid, alpha, test_data, test_classes)
            tr_refine.append(tr)
            te_refine.append(te)

        if not options.regression:
            print "Correctly classified on training data: ",tr
            print "Correctly classified on testing data:  ",te

        if options.checkpoint != None: writeCheckpoint(options.checkpoint, grid, alpha, (options.adapt_start + adaptStep))
        if options.verbose: print formTxt([te_refine[-1]], [tr_refine[-1]], [num_refine[-1]], False)

        #increment adaptive step
        adaptStep += 1
#        if(options.adaptive >= 0 and adaptStep <= options.adaptive) \
#            or (options.epochs_limit > 0 and options.epochs_limit > getEpochsErrorIncreasing(te_refine)) \
#            or (options.regression and options.mse_limit > 0 and options.mse_limit < te_refine[-1]) \
#            or (options.grid_limit > 0 and options.grid_limit > grid.getSize()) \
#            :
#        print "-----------------------------------------------------------------"
#        print options.adaptive < 0, adaptStep <= options.adaptive
#        print options.epochs_limit <= 0, options.epochs_limit > getEpochsErrorIncreasing(te_refine)
#        print options.mse_limit <= 0, options.regression, options.mse_limit < te_refine[-1]
#        print options.grid_limit <= 0 , options.grid_limit >= grid.getSize()
#        print ((options.adaptive < 0 or adaptStep <= options.adaptive) \
#            and (options.epochs_limit <= 0 or options.epochs_limit > getEpochsErrorIncreasing(te_refine)) \
#            and (options.mse_limit <= 0 or (options.regression and options.mse_limit < te_refine[-1])) \
#            and (options.grid_limit <= 0 or options.grid_limit >= grid.getSize()))
        if((options.adaptive < 0 or adaptStep <= options.adaptive) \
            and (options.epochs_limit <= 0 or options.epochs_limit > getEpochsErrorIncreasing(te_refine)) \
            and (options.mse_limit <= 0 or (options.regression and options.mse_limit < te_refine[-1])) \
            and (options.grid_limit <= 0 or options.grid_limit >= grid.getSize())) \
            :
            print("refining grid")
            if options.regression:
                grid.getGenerator().refine(SurplusRefinementFunctor(errors, getNumOfPoints(options, grid), options.adapt_threshold))
            else:
                grid.getGenerator().refine(SurplusRefinementFunctor(alpha, getNumOfPoints(options, grid), options.adapt_threshold))

            if(options.verbose): print("Number of points: %d" %(grid.getSize(),))

        #Break condition for while-loop:
        #number of adaptive steps is achieved or MSE of test data increase last 20 epochs or last MSE less then given boundary
        else:
            break

    #--end of while loop

    if options.outfile:
        writeAlphaARFF(options.outfile, alpha)
    if options.gridfile:
        writeGrid(options.gridfile, grid)

    if options.stats != None:
        txt = formTxt(te_refine, tr_refine, num_refine)
        writeStats(options.stats, txt)
        if options.verbose: print txt

    if(options.gnuplot != None):
        if(dim != 1 and dim != 2):
            print("Wrong dimension for gnuplot-Output!")
        else:
            if options.gnuplotdata:
                writeGnuplot(options.gnuplot, grid, alpha, options.res, data=training, fvals=y)
            else:
                writeGnuplot(options.gnuplot, grid, alpha, options.res)

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
                return i-1
    return length-1

##
# returns txt variable for stats
#
def formTxt(te_refine, tr_refine, num_refine, withHeader = True):
    txt = ""
    if withHeader:
        txt = "%d %-10g %2d" % (options.level, options.regparam, options.adapt_points)
    for i in xrange(len(tr_refine)):
        txt = txt + ", %d %.10f %.10f" % (num_refine[i], tr_refine[i], te_refine[i])
    return txt + "\n"

##
# returns txt variable for stats with validation set
#
def formTxtVal(te_refine, tr_refine, val_refine, num_points, withHeader = True):
    txt = ""
    if withHeader:
        txt = "%d %-10g %2d" % (options.level, options.regparam, options.adap_points)
    for i in xrange(len(tr_refine)):
        txt = txt + ", %f %.10f %.10f %.10f" % (num_points[i], tr_refine[i], val_refine[i], te_refine[i])
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
## Learn a data-set with a stratified n-fold.
# Unlike the method doFoldStratified(), here all folds will run. One can use
# doFoldStratified() with options.onlyfoldnum == -1 to achieve the same
# behaviour.
def doFoldr():
    data = openFile(options.data[0])
    (dvec, cvec) = split_n_folds_stratified(data, options.f_level, options.seed)

    if options.regression:
        performFoldRegression(dvec, cvec)
    else:
        performFold(dvec, cvec)

#-------------------------------------------------------------------------------
## Learn a dataset with a stratified n-fold.
# Unlike the method doFoldr(), doFoldStratified() is able to learn from the
# selected fold only defined in options.onlyfoldnum.
def doFoldStratified():
    data = openFile(options.data[0])
    (dvec, cvec) = split_n_folds_stratified(data, options.f_level, options.seed)

    if options.regression:
        raise Exception("Not implemented!")
    else:
        # run all folds?
        if options.onlyfoldnum == -1:
            print "Running all folds"
            for i in range(options.f_level):
                performFoldNew(dvec, cvec, i)
        # run only single one of the folds
        else:
            performFoldNew(dvec, cvec, options.onlyfoldnum)


#-------------------------------------------------------------------------------
## Learn a dataset with a n-fold from a set of files.
# In this case folding level is determined by the number of files.
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

        refinealpha = DataVector(grid.getSize())
        refinealpha.setAll(0.0)

        alpha = DataVector(refinealpha)

        for foldSetNumber in xrange(options.f_level):
#            alpha.setAll(0.0)
            training,classes = assembleTrainingVector(dvec,cvec,foldSetNumber)

            m = Matrix(grid, training, options.regparam, options.CMode, options.Hk)
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

        num_points.append(grid.getSize())
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
            grid.getGenerator().refine(SurplusRefinementFunctor(getNumOfPoints(options, grid)))

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

    # construct/split DataVectors
    training,classes = assembleTrainingVector(dvec, cvec, ifold)
    data_tr,data_val,class_tr,class_val = split_DataVectors_by_proportion_stratified(training, classes, 0.66)

    # loop until number of refinements reached
    for adaptStep in xrange(options.adaptive + 1):
        if options.verbose: print "Step %d" % (adaptStep)

        # construct and solve CG
        alpha = DataVector(grid.getSize())
        alpha.setAll(0.0)
        m = Matrix(grid, data_tr, options.regparam, options.CMode, options.Hk)
        b = m.generateb(class_tr)
        res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix,
                     options.reuse, options.verbose, max_threshold=options.max_r)
        if options.verbose: print res

        # training and validation accuracy
        tr = testVectorFast(grid, alpha, data_tr, class_tr)
        val = testVectorFast(grid, alpha, data_val, class_val)

        # compute accuracy on test set
        # Therefore construct and solve CG again for whole training data and evaluate on test data
        m = Matrix(grid, training, options.regparam, options.CMode, options.Hk)
        b = m.generateb(classes)
        res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix,
                     True, options.verbose, max_threshold=options.max_r)
        te = testVectorFast(grid, alpha, dvec[ifold], cvec[ifold])

        num_points.append(grid.getSize())
        tr_refine.append(tr)
        val_refine.append(val)
        te_refine.append(te)

        if options.verbose:
            print "num_points:", grid.getSize()
            print "training:  ", tr
            print "validating:", val
            print "testing:   ", te

        # write checkpoint
        if options.checkpoint != None:
            writeCheckpoint(options.checkpoint, grid, alpha, options.adapt_start + adaptStep, fold=ifold)

        # refine
        if(adaptStep < options.adaptive):
            if options.verbose: print "refining"
            grid.getGenerator().refine(SurplusRefinementFunctor(alpha, getNumOfPoints(options, grid)))

    # statistics
    txt = formTxtVal(te_refine, tr_refine, val_refine, num_points)
    print txt
    if options.stats != None:
        writeStats(options.stats, txt)

    return

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def performFoldRegression(dvec,cvec):
    """Perform n-fold cross-validation.
    @param dvec contains n DataMatrices for the single folds;
    @param cvec contains n DataVectors with function values for the single folds"""

    grid = constructGrid(dvec[0].getNcols())

    num_points = []
    tr_refine = []
    te_refine = []
    tr_meanSqrError = []
    te_meanSqrError = []

    for adaptStep in xrange(options.adaptive + 1):
        meanSqrErrorsTraining = []
        meanSqrErrorsTesting = []

        refineerrors = DataVector(grid.getSize())
        refineerrors.setAll(0.0)

        alpha = DataVector(grid.getSize())

        for foldSetNumber in xrange(options.f_level):
#            alpha.setAll(0.0)
            training,classes = assembleTrainingVector(dvec,cvec,foldSetNumber)

            m = Matrix(grid, training, options.regparam, options.CMode, options.Hk)
            b = m.generateb(classes)

            res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, options.reuse, options.verbose, max_threshold=options.max_r)
            print res

            # calculate squared error per basis function
            temp = DataVector(len(classes))
            m.B.mult(alpha, temp)
            temp.sub(classes)
            temp.sqr()
            # MSE for training set
            tr = temp.sum() / len(temp)
            meanSqrErrorsTraining.append(tr)
            errors = DataVector(len(alpha))
            m.B.multTranspose(temp, errors)

            # compute MSE for test set
            te = testVectorFastMSE(grid,alpha,dvec[foldSetNumber],cvec[foldSetNumber])
            meanSqrErrorsTesting.append(te)

            refineerrors.add(errors)

            if options.verbose:
                print "Fold-%d MSE (te, tr):" % (foldSetNumber), te, tr

        trSqrError = sum(meanSqrErrorsTraining)/options.f_level
        trVar = sum(map(lambda x: (x-trSqrError)**2, meanSqrErrorsTraining))/(options.f_level-1)
        teSqrError = sum(meanSqrErrorsTesting)/options.f_level
        teVar = sum(map(lambda x: (x-teSqrError)**2, meanSqrErrorsTesting))/(options.f_level-1)

        if options.verbose:
            print "testing:  ", teSqrError, teVar
            print "training: ", trSqrError, trVar

        num_points.append(grid.getSize())
        tr_meanSqrError.append(trSqrError)
        te_meanSqrError.append(teSqrError)

        refineerrors.mult(1.0/options.f_level)
        if options.checkpoint != None: writeCheckpoint(options.checkpoint, grid, refineerrors)

        if(adaptStep < options.adaptive):
            print "refine"
            grid.getGenerator().refine(SurplusRefinementFunctor(refineerrors, getNumOfPoints(options, grid)))

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
        size = size + dataset.getNrows()

    size = size - dvecs[omit].getNrows()

    training = DataMatrix(size, dvecs[0].getNcols())
    classes = DataVector(size)

    i=0
    cv = DataVector(dvecs[0].getNcols())
    for vec in dvecs:
        if vec == dvecs[omit]:
            continue
        for x in xrange(vec.getNrows()):
            vec.getRow(x, cv)
            training.setRow(i, cv)
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
# @todo remove!
# @param data a list of lists that contains the points a the training data set, coordinate-wise
# @return a instance of a DataVector that stores the training data
def buildTrainingVector(data):
    return data["data"]
#    dim = len(data["data"])
#    training = DataVector(len(data["data"][0]), dim)
#
#    # i iterates over the data points, d over the dimension of one data point
#    for i in xrange(len(data["data"][0])):
#        for d in xrange(dim):
#            training[i*dim + d] = data["data"][d][i]
#
#    return training

#-------------------------------------------------------------------------------
## Computes the classification accuracy on some test data.
# Tests on the classes {+1, -1}, cut-off at 0.
# @param grid the sparse grid
# @param alpha DataVector of surplusses
# @param test a DataVector containing a dataset of points to test on
# @param classes DataVector of correct class values
# @return classification accuracy
def testVector(grid,alpha,test,classes):
    p = DataVector(test.getNcols())
    correct = 0
    for i in xrange(test.getNrows()):
        test.getRow(i,p)
        val = createOperationEval(grid).eval(alpha,p)
        if (val < 0 and classes[i] < 0 ) or (val > 0 and classes[i] > 0 ):
            correct = correct + 1

    return float(correct)/test.getNrows()

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
    return createOperationTest(grid).test(alpha, test, classes)/float(test.getNrows())

#-------------------------------------------------------------------------------
## Computes the MSE on some test data.
#
# Tests on the classes {+1, -1}, cut-off at 0. testVectorFast uses an OpenMP enabled c++ routine for testing
# @param grid the sparse grid
# @param alpha DataVector of surplusses
# @param test a DataVector containing a dataset of points to test on
# @param vals DataVector of correct function values
# @return classification accuracy
def testVectorFastMSE(grid, alpha, test, vals):
    return createOperationTest(grid).testMSE(alpha, test, vals)


#-------------------------------------------------------------------------------
## Computes the classification accuracy on some test data.
#
# Tests on the classes {+1, -1}, cut-off at 0. testVectorFast uses an OpenMP enabled c++ routine for testing
# and displays TP TN FP FN
# @param grid the sparse grid
# @param alpha DataVector of surplusses
# @param test a DataVector containing a dataset of points to test on
# @param classes DataVector of correct class values
# @return classification accuracy
def testVectorFastWithCharacteristicNumbers(grid, alpha, test, classes):
    charaNum = DataVector(4)
    acc = 0.0
    acc = createOperationTest(grid).testWithCharacteristicNumber(alpha, test, classes, charaNum)
    print "TP: " + str(charaNum[0])
    print "TN: " + str(charaNum[1])
    print "FP: " + str(charaNum[2])
    print "FN: " + str(charaNum[3]) + " \n"
    return acc/float(test.getNrows())

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def testVectorValues(grid,alpha,test,classes,evalValues):
    p = DataVector(1,test.getDim())
    correct = 0
    for i in xrange(test.getNrows()):
        test.getRow(i,p)
        val = grid.EvaluatePoint(p, alpha)
        evalValues.append(val)
        if (val < 0 and classes[i] < 0 ) or (val > 0 and classes[i] > 0 ):
            correct = correct + 1

    return float(correct)/test.getNrows()

#-------------------------------------------------------------------------------
## Tests the classifier with some test data
#
# Tests on the classes {+1,-1}, during this tests the charateristic numbers of the classifier are computed
# : TP TN FP FN
# @param grid the sparse grid
# @param alpha DataVector of surplusses
# @param test a DataVector containing a dataset of points to test on
# @param classes DataVector of correct class values
# @param evalValues reference to an array into which the function's evaluations are stroed
def testValuesWithCharacteristicNumbers(grid,alpha,test,classes,evalValues):
    p = DataVector(1,test.getDim())
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    for i in xrange(test.getNrows()):
        test.getRow(i,p)
        val = createOperationEval(grid).eval(alpha, p)
        evalValues.append(val)
        if (val > 0 and classes[i] > 0 ):
            TP = TP + 1

        if (val < 0 and classes[i] < 0 ):
            TN = TN + 1

        if (val > 0 and classes[i] < 0 ):
            FP = FP + 1

        if (val < 0 and classes[i] > 0 ):
            FN = FN + 1

    print "TP: " + str(TP)
    print "TN: " + str(TN)
    print "FP: " + str(FP)
    print "FN: " + str(FN) + " \n"

#-------------------------------------------------------------------------------
## builds a vector that contains the class information
#
# @param data class information
# @return DataVector that contains the class information, length of the vector is identical to the length of the data vector
def buildYVector(data):
    return data["classes"]

#    y = DataVector(len(data["classes"]))
#    for i in xrange(len(data["classes"])):
#        y[i] = data["classes"][i]
#
#    return y



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
    parser.add_option("--adapt_rate", action="store", type="float", dest="adapt_rate", metavar="NUM", help="Percentage of points from all refinable points in one refinement iteration")
    parser.add_option("--adapt_start", action="store", type="int", default=0, dest="adapt_start", metavar="NUM", help="The index of adapt step to begin with")
    parser.add_option("--adapt_threshold", action="store", type="float", default=0.0, dest="adapt_threshold", metavar="NUM", help="The threshold, an error or alpha has to be greater than in order to be reined.")
    parser.add_option("-m", "--mode", action="store", type="string", default="apply", dest="mode", help="Specifies the action to do. Get help for the mode please type --mode help.")
    parser.add_option("-C", "--CMode", action="store", type="string", default="laplace", dest="CMode", help="Specifies the action to do.")
    parser.add_option("-f", "--foldlevel", action="store", type="int",default=10, metavar="LEVEL", dest="f_level", help="If a fold mode is selected, this specifies the number of sets generated")
    parser.add_option("--onlyfoldnum", action="store", type="int", default=-1, metavar="I", dest="onlyfoldnum", help="Run only fold I in n-fold cross-validation. Default: run all")
    parser.add_option("-L", "--lambda", action="store", type="float",default=0.000001, metavar="LAMBDA", dest="regparam", help="Lambda")
    parser.add_option("-i", "--imax", action="store", type="int",default=500, metavar="MAX", dest="imax", help="Max number of iterations")
    parser.add_option("-r", "--accuracy", action="store", type="float",default=0.0001, metavar="ACCURACY", dest="r", help="Specifies the accuracy of the CG-Iteration")
    parser.add_option("--max_accuracy", action="store", type="float", default=None, metavar="ACCURACY", dest="max_r", help="If the norm of the residuum falls below ACCURACY, stop the CG iterations")
    parser.add_option("-d", "--data", action="append", type="string", dest="data", help="Filename for the Datafile.")
    parser.add_option("-t", "--test", action="store", type="string", dest="test", help="File containing the testdata")
    parser.add_option("--val_proportion", action="store", type="string", dest="val_proportion", metavar="p", default=None,
                      help="Proportion (0<=p<=1) of training data to take as validation data (if applicable)")
    parser.add_option("-A", "--alpha", action="store", type="string", dest="alpha", help="Filename for a file containing an alpha-Vector")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="Filename where the calculated alphas are stored")
    parser.add_option("--gridfile", action="store", type="string", dest="gridfile", help="Filename where the resulting grid is stored")
    parser.add_option("-g", "--gnuplot", action="store", type="string", dest="gnuplot", help="In 2D case, the generated can be stored in a gnuplot readable format.")
    parser.add_option("--gnuplotdata", action="store_true", dest="gnuplotdata", default=False, help="In 2D case, the generated can be stored in a gnuplot readable format.")
    parser.add_option("-R", "--resolution", action="store", type="int",default=50, metavar="RESOLUTION", dest="res", help="Specifies the resolution of the gnuplotfile")
    parser.add_option("-s", "--stats", action="store", type="string", dest="stats", help="In this file the statistics from the test are stored")
    parser.add_option("-p", "--polynom", action="store", type="int", default=0, dest="polynom", help="Sets the maximum degree for high order basis functions. Set to 2 or larger to activate. Works only with 'identity' and 'fold'-modes.")
    parser.add_option("-b", "--border", action="store_true", default=False, dest="border", help="Enables special border base functions")
    parser.add_option("--boundary", action="store", type="int", default=False, dest="boundary", help="Use basis functions on boundary (trapezoid boundary==1, boundary==2)")
    parser.add_option("--trapezoid-boundary", action="store_true", default=False, dest="trapezoidboundary", help="Enables boundary functions that have a point on the boundary for every inner point (Trapezoid)")
    parser.add_option("--complete-boundary", action="store_true", default=False, dest="completeboundary", help="Enables boundary functions that have more points on the boundary than inner points")
    parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose", help="Provides extra output")
    parser.add_option("--normfile", action="store", type="string", dest="normfile", metavar="FILE", help="For all modes that read data via stdin. Normalizes data according to boundaries in FILE")
    parser.add_option("--reuse", action="store_true", default=False, dest="reuse", help="Reuse alpha-values for CG")
    parser.add_option("--seed", action="store", type="float", dest="seed", help="Random seed used for initializing")
    parser.add_option("--regression", action="store_true", default=False, dest="regression", help="Use regression approach.")
    parser.add_option("--checkpoint", action="store", type="string", dest="checkpoint", help="Filename for checkpointing. For fold? and test. No file extension.")
    parser.add_option("--grid", action="store", type="string", dest="grid", help="Filename for Grid-resume. For fold? and test. Full filename.")
# @todo (khakhutv) maybe a better name for parameter
    parser.add_option("--epochs_limit", action="store", type="int", default="0", dest="epochs_limit", help="Number of refinement iterations (epochs), MSE of test data have to increase, before refinement will stop.")
    parser.add_option("--mse_limit", action="store", type="float", default="0.0", dest="mse_limit", help="If MSE of test data fall below this limit, refinement will stop.")
    parser.add_option("--grid_limit", action="store", type="int", default="0", dest="grid_limit", help="If the number of points on grid exceed grid_limit, refinement will stop.")
    parser.add_option("--Hk", action="store", type="float", default="1.0", dest="Hk", help="Parameter k for regularization with H^k norm. For certain CModes.")

    parser.add_option("--function_type", action="store", type="choice", default=None, dest="function_type", choices=['modWavelet'],
                      help="Choose type for non-standard basis functions")

    # parse options
    (options,args)=parser.parse_args()

    # check some options
    CMode = options.CMode.lower()

    #this incrementation is unobvious, the default value of options.adaptive
    #options.adaptive = options.adaptive + 1

    # check C-mode
    if CMode == "help":
        print "The following C-modes are available:"
        for m in CModes.keys():
            print "%15s: %s" % (m, CModes[m])
        sys.exit(0)
    elif CMode not in CModes.keys():
        print CMode, CModes.keys()
        print("Wrong C-mode! Please refer to '-C help' for further information.")
        sys.exit(1)

    if options.polynom > 1 and CMode != "identity":
        print("Wrong C-mode selected for high-order grids.")
        sys.exit(1)

#    if options.polynom > 1 and options.border:
#        print("Special border bases do not work for High-Order Grids")
#        sys.exit(1)

    # check further parameters:
    if options.onlyfoldnum <> -1 and (not options.onlyfoldnum in range(options.f_level)):
        parser.error("--onlyfoldnum: Not in range 0,...,--foldlevel -1")

    # adapt_rate has to be <1
    if options.adapt_rate != None and (options.adapt_rate <= 0 or options.adapt_rate >= 1):
        parser.error("--adapt_rate has to be between 0 and 1")



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
        'foldstratified'    : {'help': "learn a dataset with a stratified n-fold using a validation set",
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
