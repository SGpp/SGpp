from bin.data.DataContainer import DataContainer
from bin.uq.analysis import KnowledgeTypes
from bin.uq.learner.Learner import Learner, LearnerEvents
from bin.uq.sampler import Samples, SampleType
from pysgpp import DataVector, SurplusRefinementFunctor, DataMatrix

import numpy as np


class Regressor(Learner):
    """
    Subclass of Learner, responsible for regression.
    The methods specific for regression are implemented here.
    @todo (khakhutv) implement a test case
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
        refinableNum = self.grid.createGridGenerator().getNumberOfRefinablePoints()
        pointsNum = self.getNumOfPointsToRefine(refinableNum)
        # @todo (khakhutv) (low) develop a way to simplify interfaces and
        # use different functors
        functor = SurplusRefinementFunctor(self.errors, pointsNum,
                                           self.getAdaptThreshold())
        self.grid.createGridGenerator().refine(functor)

#     ## Creates Regressor object without any parameter set
#     # @param jsonObject: A json object.
#     # @return: Regressor object without any parameter set
#    @classmethod
#    def fromJson(cls, jsonObject):
#        return Regressor()
