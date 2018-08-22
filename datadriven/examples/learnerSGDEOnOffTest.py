import numpy as np

import pysgpp as sg
from pysgpp.extensions.datadriven.learner import Types

###
# This example shows how to perform offline/online-classification using sparse
# grid density estimation and matrix decomposition methods. It creates an
# instance of LearnerSGDEOnOff and runs the function train() where the
# main functionality is implemented.
#
# Currently, only binary classification with class labels -1 and 1 is possible.
#
# The example provides the option to execute several runs over differently
# ordered data and perform a 5-fold cross-validation within each run.
# Therefore,
# already randomly ordered and partitioned data is required.
# Average results from several runs might be more reliable in an
# online-learning
# scenario, because the ordering of the data points seen by the learner
# can affect the result.
##

def main():
    
    ###
    # Specify the number of runs to perform.
    # If only one specific example should be executed, set
    # totalSets=1
    ##
    totalSets = 1
    totalFolds = 1   # set to 5 to perform 5-fold cross-validation
    avgError = 0.0
    avgErrorFolds = 0.0
    
    for numSets in range(totalSets):
        
        ###
        # Vector to compute average classification error throughout the
        # learning process. The length of the vector determines the total
        # number of error observations
        ##
        avgErrorsFolds = sg.DataVector(51, 0.0)
    
        for numFolds in range(totalSets):
            
            ###
            # Get the training, test and validation data
            ##
            arff_provider = sg.ArffFileSampleProvider()

            filename = "../tests/data/ripleyGarcke.train.arff"
            print "# loading training samples: ", filename
            arff_provider.readFile(filename)
            trainDataset = arff_provider.getAllSamples()

            filename = "../tests/data/ripleyGarcke.test.arff"
            print "# loading test samples: ", filename
            arff_provider.readFile(filename)
            testDataset = arff_provider.getAllSamples()

            ###
            # Specify the number of classes and the corresponding class labels
            ##
            classNum = 2
            classLabels = sg.DataVector(classNum)
            classLabels[0] = -1
            classLabels[1] = 1

            ###
            # Configure grid
            ##
            print "# create grid config... ",
            gridConfig = sg.RegularGridConfiguration()
            gridConfig.dim_ = trainDataset.getDimension()
            gridConfig.level_ = 3
            gridConfig.type_ = sg.GridType_Linear
            print "Done"

            ###
            # Configure regularization
            ##
            print "# create regularization config... ",
            regularizationConfig = sg.RegularizationConfiguration()
            regularizationConfig.regType_ = sg.RegularizationType_Identity
            print "Done"

            ###
            # Select the desired decomposition type for the offline step
            # NOTE: Refinement/Coarsening only possible for Cholesky decomposition
            ##
            dt = sg.DBMatDecompostionType_DenseIchol
            decompType = "Incomplete Cholesky decomposition on Dense Matrix"
            print "Decomposition type: ", decompType

            ###
            # Configure adaptive refinement (if Cholesky is chosen!). As refinement
            # monitor the periodic monitor or the convergence monitor can be chosen.
            # Possible refinement indicators are surplus refinement, data-based
            # refinement, zero-crossing-based refinement
            ##
            print "# create adaptive refinement config"
            # select periodic monitor - perform refinement in fixed intervals
            refMonitor = "periodic"
            # select convergence monitor - perform refinements if algorithm has
            # converged (convergence measured with respect to changes of the
            # classification accuracy)
            #refMonitor = "convergence"

            refPeriod = 40   # refinement interval
            accDeclineThreshold = 0.001   # convergence threshold
            accDeclineBufferSize = 140   # no. of accuracy measurements which
                                         # are considered to be performed
            minRefInterval = 10   # minimum no. of iterations before
                                  # next refinement is allowed to be performed
            print "Refinement monitor: ", refMonitor

            # select zero-crossing-based refinement
            refType = "zero"
            # select surplus refinement
            #refType = "surplus"
            # select data-based refinement
            #refType = "data"
            print "Refinement type: ", refType

            ###
            # Specify number of refinement steps and the maximum number of grid
            # points to refine each step
            ##
            adaptConfig = sg.AdaptivityConfiguration()
            adaptConfig.numRefinements_ = 2
            adaptConfig.noPoints_ = 7
            adaptConfig.threshold_ = 0.0   # only required for surplus refinement!

            lambda_ = 0.01   # initial regularization parameter
            beta = 0.0   # initial weigting factor

            dconf = sg.DBMatDensityConfiguration(gridConfig, adaptConfig,
                                                 regularizationConfig.regType_, 
                                                 lambda_, dt)

            usePrior = False   # specify if prior should be used to prdict classes

            dconf.icholParameters.sweepsDecompose = 2
            dconf.icholParameters.sweepsRefine = 2

            ###
            # Create the learner
            ##
            print "# create learner... ",
            learner = sg.LearnerSGDEOnOff(dconf, trainDataset, testDataset, None,
                                          classLabels, classNum, usePrior, beta, lambda_)
            print "Done"

            ###
            # Configure cross-validation
            # Set enableCV=True to perform cross-validation during the learning process
            ##
            enableCV = False
            if enableCV:
                print "# create cross-validation config"

            # set cross-validation config, if cross-validation is enabled
            nextCVStep = 50
            cvLambdaStart = 1e-1
            cvLambdaEnd = 1e-10
            cvLambdaSteps = 10
            cvLogScale = True
            cvTestData = testDataset.getData()
            learner.setCrossValidationParameters(cvLambdaSteps, cvLambdaStart, cvLambdaEnd,
                                                 cvTestData, None, cvLogScale)

            ###
            # Learn the data
            ##
            batchSize = 1   # specify batch size (set to 1 for processing only a single
                            # data point each iteration)
            maxDataPasses = 2   # specify maximum no. of passes over training data set

            print "# start to train the learner..."
            learner.train(batchSize, maxDataPasses, refType, refMonitor, refPeriod, 
                          accDeclineThreshold, accDeclineBufferSize, minRefInterval,
                          enableCV, nextCVStep)

            ###
            # Compute accuracy on test data
            ##
            acc = learner.getAccuracy()
            print "# accuracy (test data): ", acc

            # store results (classified data, grids, density function)
            #learner.storeResults()

            tmp = sg.DataVector(0)
            avgErrorFolds += 1.0 - acc
            learner.getAvgErrors(tmp)
            avgErrorsFolds.add(tmp)
          
        ###
        # Average accuracy on test data regarding 5-fold cross-validation
        ##
        avgErrorFolds = avgErrorFolds / float(totalFolds)
        
        if totalSets > 1 and totalFolds > 1:
            s = "Average accuracy on test data (set {}):".format(numSets + 1)
            print s, (1.0 - avgErrorFolds)
        
        avgError += avgErrorFolds
        avgErrorFolds = 0.0   # reset avgErrorFolds!
        avgErrorsFolds.mult(1.0 / float(totalFolds))


if __name__ == '__main__':
    main()


