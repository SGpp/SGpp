# This file is part of SGClass, a program package making use of spatially adaptive sparse grids to solve numerical problems.
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
# @todo Join the different modes so that not 100 different versions (checkpointing, ...) exist!


from optparse import OptionParser
import sys
from tools import *
from pysgpp import *
from painlesscg import cg,sd,cg_new
from math import sqrt
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
        print ("More on the usage of %s with --help" % (sys.argv[0]))
        sys.exit(1)

    # execute action
    modes[mode]['action']()
    


#-------------------------------------------------------------------------------
## Opens and read the data of an ARFF file
# Reads in an ARFF file
#
# @param filename of the file
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
                #grid = SpGridHighOrder(dim,options.level,options.polynom)
                grid = Grid.createPolyGrid(dim, options.polynom)
            else:
                if options.verbose:
                    print "LinearGrid, l=%s" % (options.level)
                #grid = SpGridLinear(dim,options.level)
                grid = Grid.createLinearGrid(dim)
        generator  = grid.createGridGenerator()
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
    
#    # construct corresponding regular grid by dim and level
#    grid = SpGridLinear(dim, options.level)
    # construct corresponding grid
    grid = constructGrid(dim)
    if(alpha.getSize() != grid.getStorage().size()):
        print "Error: Inconsistent Grid and Alpha-Vector"
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
            
    # get filename for output file
    if(options.outfile != None):
        data["filename"] = options.outfile
    else:
        data["filename"] = data["filename"] + ".out.arff"
    
    # write data with function values to ARFF-file
    data["classes"] = classes
    writeDataARFF([data])


#-------------------------------------------------------------------------------
## Evaluate a sparse grid function given by grid and alphas at data points.
# Evaluates a sparse grid function, given by a regular sparse grid (dim, level)
# and an alpha-Vector (ARFF-file) at a set of data points (ARFF-file).
# Outputs function values to ARFF-file
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
    #print "Classifier dim: %d" % dim
#    print "Classifier q.getDim(): %d" % q.getDim()
    
    #q = [1.0,1.0,1.0,1.0,1.0]
    classes = []
    # if test data contains function values, additionally compute L2-norm of error
    if data.has_key("classes"):
        compute_accuracy = True
        err = 0
    else:
        compute_accuracy = False
    # traverse Data
    for i in xrange(numData):
        x.getRow(i,q)
#        print "q dim in DataVector: %d" % q.getDim()
#        print "q:"
#        for j in xrange(q.getSize()):
#            print q[j]
        #TODO: replace the parameters everywhere
        val = grid.createOperationEval().eval(alpha,q)
        #print "val %f" % val
        classes.append(val)
    # output accuracy:
    if compute_accuracy:
        error = DataVector(len(classes))
        for i in range(len(classes)):
            error[i] = classes[i]
#            print str(error[i]) + " - " + str(y[i])
        
        error.sub(y) # error vector
        error.sqr()  # entries squared
        # output some statistics
        err_min = error.min(1)
        err_max = error.max(1)
        print "(Min,Max) error: (%f,%f)" % (sqrt(err_min), sqrt(err_max))
        # output accuracy
        err = error.sum()
        print "L2-norm of error on data: %f" % (sqrt(err))
        print "MSE: %f" %(err / error.getSize())

    # get filename for output file
    if(options.outfile != None):
        data["filename"] = options.outfile
    else:
        data["filename"] = data["filename"] + ".out.arff"
    # write data with function values to ARFF-file
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

    for adaptStep in xrange(options.adaptive):
        alpha = DataVector(grid.getStorage().size())
        alpha.setAll(0.0)

        m = Matrix(grid, training, options.regparam, options.zeh)
        b = m.generateb(classes)
            
        res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, False, options.verbose)
        print res

        if options.regression:
            error = DataVector(len(classes))
            m.B.multTranspose(alpha, m.x, error)
            error.sub(classes)# error vector
            error.sqr()       # entries squared
            # output some statistics
            err_min = error.min(1)
            err_max = error.max(1)
            print "(Min,Max) error: (%f,%f)" % (sqrt(err_min), sqrt(err_max))
            # output accuracy
            err = error.sum()
            print "L2-norm of error on data: %f" % (sqrt(err))
            print "MSE: ",(err / error.getSize())

            # calculate error per basis function
            errors = DataVector(alpha.getSize())
            m.B.mult(error, m.x, errors)

        if options.checkpoint != None:
            txt = ""
            writeCheckpoint(options.checkpoint, grid, alpha, txt, adaptStep)
        
        
        if(adaptStep + 1 < options.adaptive):
            #if(options.verbose):
            print("refining grid")
            if options.regression:
                grid.createGridGenerator().refine(SurplusRefinementFunctor(errors))
            else:
                grid.createGridGenerator().refine(SurplusRefinementFunctor(alpha))
    return alpha

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

    for adaptStep in xrange(options.adaptive):
        m = Matrix(grid, training, options.regparam, options.zeh)
        b = m.generateb(y)
        
        alpha = DataVector(grid.getStorage().size())
        alpha.setAll(0.0)
        res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, False, options.verbose)
        print res

        tr = testVectorFast(grid, alpha, training, y)
        
#        if options.verbose:
#            teValues = []
#            te = testVectorValues(grid, alpha, test_data, test_classes, teValues)
#        else:
#            te = testVector(grid, alpha, test_data, test_classes)
        
        te = testVectorFast(grid, alpha, test_data, test_classes)


        te_refine.append(te)
        tr_refine.append(tr)
        num_refine.append(grid.getStorage().size())

        if options.verbose:
            print "Correct classified on training data: ",tr
            print "Correct classified on testing data:  ",te

        if options.checkpoint != None:
            txt = "%f, %-10g, %f" % (options.level, options.regparam, options.adaptive)
            for i in xrange(len(tr_refine)):
                txt = txt + ", %f, %.10f, %.10f" % (num_refine[i], tr_refine[i], te_refine[i])
            
            writeCheckpoint(options.checkpoint, grid, alpha, txt, adaptStep)

        if(adaptStep + 1 < options.adaptive):
            if(options.verbose):
                print("refining grid")
            grid.createGridGenerator().refine(SurplusRefinementFunctor(alpha))
            if(options.verbose):
                print("Number of points: %d" %(grid.getStorage().size(),))

    if options.stats != None:
        txt = "%f, %-10g, %f" % (options.level, options.regparam, options.adaptive)
        for i in xrange(len(tr_refine)):
            txt = txt + ", %f, %.10f, %.10f" % (num_refine[i], tr_refine[i], te_refine[i])
        if options.verbose:
            print txt
        writeLockFile(options.stats, txt+"\n")



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
    
#    if options.grid == None:
#        if options.polynom > 1:
#            grid = SpGridHighOrder(dvec[0].getDim(),options.level,options.polynom)
#        else:
#            grid = SpGridLinear(dvec[0].getDim(),options.level)
#    else:
#        grid = readGrid(options.grid)
#        
#    if options.border:
#        grid.setUseBorderFunctions(True)

    grid = constructGrid(dvec[0].getDim())
        
    num_points = []
    tr_refine = []
    te_refine = []
    
    for adaptStep in xrange(options.adaptive):
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

            res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, options.reuse, options.verbose)
            print res
            
#            tr = testVector(grid,alpha,training,classes)
#            te = testVector(grid,alpha,dvec[foldSetNumber],cvec[foldSetNumber])

            tr = testVectorFast(grid, alpha, training, classes)
            te = testVectorFast(grid, alpha, dvec[foldSetNumber], cvec[foldSetNumber])

            trainingCorrect.append(tr)
            testingCorrect.append(te)

## Verstehe ich nicht. Habs deshalb auskommentiert (Dirk)
##             if(adaptStep +1 == options.adaptive):
##                 #Letzte verfeinerung, wir sind fertig
##                 pass
##             else:
##                 refinealpha.add(alpha)
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
        
        if options.checkpoint != None:
            txt = "%f, %-10g, %f" % (options.level, options.regparam, options.adaptive)
            for i in xrange(len(tr_refine)):
                txt = txt + ", %f, %.10f, %.10f" % (num_refine[i], tr_refine[i], te_refine[i])
            
            writeCheckpoint(options.checkpoint, grid, refinealpha, txt)
        
        if options.verbose:
            print "alpha"
            print refinealpha
            print "grid"
            print grid
        if(adaptStep + 1 < options.adaptive):
            print "refine"
            #grid.refineOneGridPoint(refinealpha)
            grid.createGridGenerator().refine(SurplusRefinementFunctor(refinealpha))

    #print(tr)
    #print(te)

    if options.stats != None:
        txt = "%f, %-10g, %f" % (options.level, options.regparam, options.adaptive)
        for i in xrange(len(tr_refine)):
            txt = txt + ", %f, %.10f, %.10f" % (num_points[i], tr_refine[i], te_refine[i])
        if options.verbose:
            print txt
        writeLockFile(options.stats, txt+"\n")

    return

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
def performFoldRegression(dvec,cvec):
    
#    if options.polynom > 1:
#        grid = SpGridHighOrder(dvec[0].getDim(),options.level,options.polynom)
#    else:
#        grid = SpGridLinear(dvec[0].getDim(),options.level)
#        
#    if options.border:
#        grid.setUseBorderFunctions(True)
    
    grid = constructGrid(dvec[0].getDim())
        
    num_points = []
    tr_refine = []
    te_refine = []
    tr_meanSqrError = []
    te_meanSqrError = []
        
    for adpatStep in xrange(options.adaptive):
        trainingCorrect = []
        testingCorrect = []
        meanSqrErrorsTraining = []
        meanSqrErrorsTesting = []        

        #refinealpha = DataVector(grid.getPointsCount())
        #refinealpha.setAll(0.0)

        refineerrors = DataVector(grid.getStorage().size())
        refineerrors.setAll(0.0)

        alpha = DataVector(grid.getStorage().size())
        
        for foldSetNumber in xrange(options.f_level):
#            alpha.setAll(0.0)
            training,classes = assembleTrainingVector(dvec,cvec,foldSetNumber)
            
            m = Matrix(grid, training, options.regparam, options.zeh)
            b = m.generateb(classes)

            res = cg_new(b, alpha, options.imax, options.r, m.ApplyMatrix, options.reuse, options.verbose)
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
                #refinealpha.add(alpha)
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
        
        #refinealpha.mult(1.0/options.f_level)
        refineerrors.mult(1.0/options.f_level)
        
        if(adpatStep + 1 < options.adaptive):
            print "refine"
            #grid.refineOneGridPoint(refinealpha)
            grid.createGridGenerator().refine(SurplusRefinementFunctor(refineerrors))

    #print(tr)
    #print(te)

    if options.stats != None:
        txt = "%f, %-10g, %f" % (options.level, options.regparam, options.adaptive)
        for i in xrange(len(tr_refine)):
            txt = txt + ", %f, %.10f, %.10f" % (num_points[i], tr_meanSqrError[i], te_meanSqrError[i])
        if options.verbose:
            print txt
        writeLockFile(options.stats, txt+"\n")

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
## Computes the classification accuracy on some test data. Tests on the classes
# {+1, -1}, cut-off at 0.
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
#-------------------------------------------------------------------------------
## Computes the classification accuracy on some test data. Tests on the classes
# {+1, -1}, cut-off at 0. testVectorFast uses an OpenMP enabled c++ routine for testing
# @param grid the sparse grid
# @param alpha DataVector of surplusses
# @param test a DataVector containing a dataset of points to test on
# @param classes DataVector of correct class values
# @return classification accuracy
def testVectorFast(grid, alpha, test, classes):
    return float(grid.createOperationEval().test(alpha, test, classes))/test.getSize()

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
#-------------------------------------------------------------------------------
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
    parser.add_option("-a", "--adaptive", action="store", type="int", default="0", dest="adaptive", metavar="NUM", help="Using an adaptive Grid with NUM of refines")
    parser.add_option("-m", "--mode", action="store", type="string", default="apply", dest="mode", help="Specifies the action to do. Get help for the mode please type --mode help.")
    parser.add_option("-C", "--zeh", action="store", type="string", default="laplace", dest="zeh", help="Specifies the action to do.")
    parser.add_option("-f", "--foldlevel", action="store", type="int",default="10", metavar="LEVEL", dest="f_level", help="If a fold mode is selected, this specifies the number of sets generated")
    parser.add_option("-L", "--lambda", action="store", type="float",default="0.000001", metavar="LAMBDA", dest="regparam", help="Lambda")
    parser.add_option("-i", "--imax", action="store", type="int",default="400", metavar="MAX", dest="imax", help="Max number of iterations")
    parser.add_option("-r", "--accuracy", action="store", type="float",default="0.0001", metavar="ACCURACY", dest="r", help="Specifies the accuracy of the CG-Iteration")
    parser.add_option("-d", "--data", action="append", type="string", dest="data", help="Filename for the Datafile.")
    parser.add_option("-t", "--test", action="store", type="string", dest="test", help="File containing the testdata")
    parser.add_option("-A", "--alpha", action="store", type="string", dest="alpha", help="Filename for a file containing an alpha-Vector")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="Filename where the calculated alphas are stored")
    parser.add_option("-g", "--gnuplot", action="store", type="string", dest="gnuplot", help="In 2D case, the generated can be stored in a gnuplot readable format.")
    parser.add_option("-R", "--resolution", action="store", type="int",default="50", metavar="RESOLUTION", dest="res", help="Specifies the resolution of the gnuplotfile")
    parser.add_option("-s", "--stats", action="store", type="string", dest="stats", help="In this file the statistics from the test are stored")
    parser.add_option("-p", "--polynom", action="store", type="int", default="0", dest="polynom", help="Sets the maximum degree for high order basis functions. Set to 2 or larger to activate. Works only with 'identity' and 'fold'-modes.")
    parser.add_option("-b", "--border", action="store_true", default=False, dest="border", help="Enables special border base functions")
    parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose", help="Provides extra output")
    parser.add_option("--normfile", action="store", type="string", dest="normfile", metavar="FILE", help="For all modes that read data via stdin. Normalizes data according to boundaries in FILE")
    parser.add_option("--reuse", action="store_true", default=False, dest="reuse", help="Reuse alpha-values for CG")
    parser.add_option("--seed", action="store", type="int", dest="seed", help="Random seed used for initializing")
    parser.add_option("--regression", action="store_true", default=False, dest="regression", help="Use regression approach.")
    parser.add_option("--checkpoint", action="store", type="string", dest="checkpoint", help="Filename for checkpointing. For fold? and test. No file extension.")
    parser.add_option("--grid", action="store", type="string", dest="grid", help="Filename for Grid-resume. For fold? and test. Full filename.")
    # parse options
    (options,args)=parser.parse_args()

    # check some options
    zeh = options.zeh.lower()

    options.adaptive = options.adaptive + 1

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
        'foldr'    : {'help': "learn a dataset with a stratified n-fold",
                      'required_options': ['data', ['level', 'grid']],
                      'action': doFoldr},
        'foldf'    : {'help': "learn a dataset with a n-fold from a set of files",
                      'required_options': ['data', ['level', 'grid']],
                      'action': doFoldf}
        }

    # Execute the mode
    exec_mode(options.mode.lower())
