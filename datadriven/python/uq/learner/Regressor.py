from pysgpp.extensions.datadriven.data.DataContainer import DataContainer
from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.learner.Learner import Learner, LearnerEvents
from pysgpp.extensions.datadriven.uq.sampler import Samples, SampleType
from pysgpp import DataVector, SurplusRefinementFunctor, DataMatrix

import numpy as np


class Regressor(Learner):
    """
    Subclass of Learner, responsible for regression.
    The methods specific for regression are implemented here.
    """

    def __init__(self):
        """
        Constructor
        """
        super(self.__class__, self).__init__()
        # Errors per basis function
        self.errors = None
        # Error vector
        self.error = None

    def __getattr__(self, attr):
        """
        Overrides built-in method if method called is not a object
        method of this Descriptor, most probably it's a method of
        the learner so it tries to call the method
        from our specification
        @param attr: string method name
        @return: method call in specification
        """
        return getattr(self.specification, attr)

    # # Learn data from training data set and use validation data set to prevent overfitting
    #
    # @param dataset: DataContainer object with data sets, default value None (initialized data set used)
    # @return: DataVector of alpha
    def learnDataWithTest(self, dataset=None):
        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_STARTED)
        B = createOperationMultipleEval(self.grid,
                                        self.dataContainer.getPoints(DataContainer.TRAIN_CATEGORY))
        self.specification.setBOperator(B)

        if dataset is None:
            dataset = self.dataContainer

        # learning step
        trainSubset = dataset.getTrainDataset()
        # testpoint = data.allPoint\points
        # testvalues = data.allValues\values
        testSubset = dataset.getTestDataset()

        while True:  # repeat until policy says "stop"
            self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_STEP_STARTED)

            self.alpha = self.doLearningIteration(trainSubset)

            # calculate avg. error for training and test data and avg. for refine alpha
            self.updateResults(self.alpha, trainSubset, testSubset)

            self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_STEP_COMPLETE)

            self.iteration += 1

            if self.stopPolicy.isTrainingComplete(self):
                break

            # refine grid
            self.refineGrid()

        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_COMPLETE)
        return self.alpha

    # # Simple data learning
    #
    # @return: DataVector of alpha
    def learnData(self):
        self.notifyEventControllers(LearnerEvents.LEARNING_STARTED)
        self.specification.setBOperator(createOperationMultipleEval(self.grid,
                    self.dataContainer.getPoints(DataContainer.TRAIN_CATEGORY)))
        print self.getL()
        while True:  # repeat until policy says "stop"
            print "Learning %i/%i" % (self.iteration, self.stopPolicy.getAdaptiveIterationLimit())
            self.notifyEventControllers(LearnerEvents.LEARNING_STEP_STARTED)
            # learning step
            self.alpha = self.doLearningIteration(self.dataContainer)

            # calculate avg. error for training and test data and avg. for refine alpha
            self.updateResults(self.alpha, self.dataContainer)
            self.notifyEventControllers(LearnerEvents.LEARNING_STEP_COMPLETE)
            self.iteration += 1
            if (self.stopPolicy.isTrainingComplete(self)):
                break
            # refine grid
            self.refineGrid()

#         from pysgpp.extensions.datadriven.uq.plot import plotNodal3d
#         plotNodal3d(self.grid, self.alpha)
#         data = self.dataContainer.getPoints('train').array()
#         fig = plt.figure()
#         plt.plot(data[:, 0], data[:, 1], ' ', marker='v')
#         fig.show()
#         plt.show()

        self.notifyEventControllers(LearnerEvents.LEARNING_COMPLETE)
        return self.alpha

    # # Learn data with cross-fold validation
    #
    # @return: list of DataVector alpha in different folds
    def learnDataWithFolding(self,):
        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_FOLDING_STARTED)
        self.specification.setBOperator(createOperationMultipleEval(self.grid,
                  self.dataContainer.getPoints(DataContainer.TRAIN_CATEGORY)))
        # update folding
        self.updateFoldingPolicy()
        alphas = []
        for dataset in self.foldingPolicy:
            alphas.append(self.learnDataWithTest(dataset))

        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_FOLDING_COMPLETE)
        return alphas

    # # Perform one learning step
    #
    # @param set: DataContainer training data set
    # @return: DataVector alpha vector
    def doLearningIteration(self, set):
        # initialize values
        self.linearSystem = DMSystemMatrix(self.grid,
                                           set.getPoints(),
                                           self.specification.getCOperator(),
                                           self.specification.getL())
        size = self.grid.getSize()
        # Reuse data from old alpha vector increasing its dimension
        if self.solver.getReuse() and self.alpha is not None:
            alpha = DataVector(self.alpha)
            alpha.resize(size)
        # Use new alpha vector
        else:
            alpha = DataVector(size)
            alpha.setAll(0.0)
        b = DataVector(size)
        self.linearSystem.generateb(set.getValues(), b)
        # calculates alphas
        self.solver.solve(self.linearSystem, alpha, b, self.solver.getReuse(),
                          False, self.solver.getThreshold())

        return alpha

    def getL2NormError(self):
        """
        calculate L2-norm of error
        @return: last L2-norm of error
        """
        return np.sqrt(self.error.sum())

    def getMaxError(self):
        """
        calculate max error
        @return: max error
        """
        return np.sqrt(self.error.max())

    def getMinError(self):
        """
        calculate min error
        @return: min error
        """
        return np.sqrt(self.error.min())

    def evalError(self, data, alpha):
        """
        Evaluate regression MSE
        @param data: DataContainer data set
        @param alpha: DataVector alpha-vector
        @return: mean square error
        """
        size = data.getPoints().getNrows()
        if size == 0:
            return 0

        self.error = DataVector(size)
        self.getBOperator().mult(alpha, self.error)
        # error vector
        self.error.sub(data.getValues())
        # entries squared
        self.error.sqr()
        errorsum = self.error.sum()
        # MSE
        mse = errorsum / size

        # calculate error per basis function
        self.errors = DataVector(len(alpha))
        self.getBOperator().multTranspose(self.error, self.errors)
        self.errors.componentwise_mult(alpha)

        # calculate error per basis function
#        self.errors = DataVector(alpha.getSize())
#        self.specification.getBOperator().mult(self.error, data.getPoints(), self.errors)

        return mse

    def updateResults(self, alpha, trainSubset, testSubset=None):
        """
        Update different statistics about training progress
        @param alpha: DataVector alpha-vector
        @param trainSubset: DataContainer with training data
        @param testSubset: DataContainer with validation data
        """
        # self.knowledge.update(alpha)
        # eval Error for training data and append it to other in this iteration
        self.trainAccuracy.append(self.evalError(trainSubset, alpha))

        i = float(len(self.trainAccuracy))

        # eval error for test data and append it to other in this iteration
        if testSubset is not None:
            self.testAccuracy.append(self.evalError(testSubset, alpha))
            self.testingOverall.append(sum(self.testAccuracy) / i)

        self.trainingOverall.append(sum(self.trainAccuracy) / i)

        self.numberPoints.append(self.grid.getSize())

    def refineGrid(self):
        """
        Refines grid with the number of points as specified in corresponding
        TrainingSpecification object
        """
        self.notifyEventControllers(LearnerEvents.REFINING_GRID)
        refinableNum = self.grid.getGenerator().getNumberOfRefinablePoints()
        pointsNum = self.getNumOfPointsToRefine(refinableNum)
        functor = SurplusRefinementFunctor(self.errors, pointsNum,
                                           self.getAdaptThreshold())
        self.grid.getGenerator().refine(functor)

#     ## Creates Regressor object without any parameter set
#     # @param jsonObject: A json object.
#     # @return: Regressor object without any parameter set
#    @classmethod
#    def fromJson(cls, jsonObject):
#        return Regressor()
