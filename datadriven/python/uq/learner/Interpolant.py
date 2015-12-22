from pysgpp.extensions.datadriven.uq.operations import (hierarchize,
                               evalSGFunctionMulti,
                               hierarchizeBruteForce)
from pysgpp import (DataVector, DataMatrix)

from Learner import Learner, LearnerEvents
import numpy as np
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import checkInterpolation, dehierarchize
from pysgpp.extensions.datadriven.uq.plot import plotNodal3d, plotSGNodal3d
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotSG3d
import pdb
from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes


class Interpolant(Learner):

    def __init__(self):
        super(Interpolant, self).__init__()
        # Error per basis function
        self.errors = None
        # Error vector
        self.error = None

        # init stats
        self.trainAccuracy = []
        self.trainCount = []
        self.trainingOverall = []
        self.testAccuracy = []
        self.testCount = []
        self.testingOverall = []

    def doLearningIteration(self, points):
        """
        Interpolates the given points with the current grid
        @param points: interpolation points
        @return: Return hierarchical surpluses
        """
        gs = self.grid.getStorage()

        # assert that the number of dimensions of the data is the same
        # as the grids
        assert gs.dim() == points.getDim()

        nodalValues = DataVector(gs.size())
        nodalValues.setAll(0.0)

        # interpolation on nodal basis
        p = DataVector(gs.dim())
        cnt = 0
        for i in xrange(gs.size()):
            gp = gs.get(i)
            gp.getCoords(p)
            x = tuple(p.array())
            if x not in points:
                # # search for 2*d closest grid points
                # q = DataVector(gs.dim())
                # l = np.array([])
                # for j in xrange(gs.size()):
                #     gs.get(j).getCoords(q)
                #     q.sub(p)
                #     l = np.append(l, q.l2Norm())

                # n = min(gs.size(), gs.dim())

                # ixs = np.argsort(l)
                # # nodalValues[i] = np.mean(l[ixs[:n]])
                nodalValues[i] = 0.0
                print p, nodalValues[i]
                cnt += 1
            else:
                nodalValues[i] = float(points[x])

        if cnt > 0:
            print '%i/%i of the grid points have \
                   been set to 0' % (cnt, gs.size())
            pdb.set_trace()

        # hierarchization
        alpha = hierarchize(self.grid, nodalValues)

        # -----------------------------------------
        # check if interpolation property is given
#         fig, _ = plotNodal3d(A)
#         fig.show()
#         fig, _ = plotSGNodal3d(self.grid, alpha)
#         fig.show()
#         fig, _ = plotSG3d(self.grid, alpha)
#         fig.show()

        err, _ = checkInterpolation(self.grid, alpha, nodalValues, epsilon=1e-12)

        if len(err) > 0:
            print "interpolation property not met"
            pdb.set_trace()
        # -----------------------------------------

        return alpha

    def learnData(self, *args, **kws):
        # learning step
        self.notifyEventControllers(LearnerEvents.LEARNING_STARTED)
        # load data sets
        trainSubset = self.dataContainer.getTrainDataset()

        # learning step
        self.alpha = self.doLearningIteration(trainSubset)
        self.updateResults(self.alpha, trainSubset, *args, **kws)

        self.iteration += 1
        self.notifyEventControllers(LearnerEvents.LEARNING_COMPLETE)
        return self.alpha

    def learnDataWithTest(self, dataset=None, *args, **kws):
        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_STARTED)

        if dataset is None:
            dataset = self.dataContainer

        # load data sets
        trainSubset = dataset.getTrainDataset()
        testSubset = dataset.getTestDataset()

        # learning step
        self.alpha = self.doLearningIteration(trainSubset)
        self.updateResults(self.alpha, trainSubset, testSubset,
                           *args, **kws)
        self.iteration += 1
        self.notifyEventControllers(LearnerEvents.LEARNING_WITH_TESTING_COMPLETE)
        return self.alpha

    def updateResults(self, alpha, trainSubset, testSubset=None,
                      dtype=KnowledgeTypes.SIMPLE):
        # evaluate MSE of training data set -> should be zero
        mse = self.evalError(trainSubset, alpha)
        self.trainAccuracy.append(mse)
        self.trainCount.append(len(self.error))
        self.trainingOverall.append(float(np.mean(self.trainAccuracy)))

        # store interpolation quality
        if testSubset:
            # L2 error and MSE
            mse = self.evalError(testSubset, alpha)
            self.testAccuracy.append(mse)
            self.testCount.append(len(self.error))
            # testing overall
            self.testingOverall.append(float(np.mean(self.testAccuracy)))

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
        computes the discrete L2-error of the sparse grid interpolant
        with respect to some MC samples at given time steps
        @param data: DataContainer samples
        @param alpha: DataVector hierarchical coefficients
        @return: mean squared error
        """
        points = data.getPoints()
        values = data.getValues()
        size = points.getNrows()
        if size == 0:
            return 0

        self.error = evalSGFunctionMulti(self.grid, alpha, points)
        # compute L2 error
        self.error.sub(values)
        self.error.sqr()

        return self.error.sum() / len(self.error)
