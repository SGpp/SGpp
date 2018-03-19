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
        self.trainErrors = np.array([])
        self.testErrors = np.array([])

    def doLearningIteration(self, points):
        """
        Interpolates the given points with the current grid
        @param points: interpolation points
        @return: Return hierarchical surpluses
        """
        gs = self.grid.getStorage()

        # assert that the number of dimensions of the data is the same
        # as the grids
        assert gs.getDimension() == points.getDim()

        nodalValues = np.ndarray(gs.getSize())

        # interpolation on nodal basis
        p = DataVector(gs.getDimension())
        cnt = 0
        for i in xrange(gs.getSize()):
            gp = gs.getPoint(i)
            gs.getCoordinates(gp, p)
            x = tuple(p.array())
            if x not in points:
                # # search for 2*d closest grid points
                # q = DataVector(gs.getDimension())
                # l = np.array([])
                # for j in xrange(gs.getSize()):
                #     gs.getCoordinates(gs.getPoint(j), q)
                #     q.sub(p)
                #     l = np.append(l, q.l2Norm())

                # n = min(gs.getSize(), gs.getDimension())

                # ixs = np.argsort(l)
                # # nodalValues[i] = np.mean(l[ixs[:n]])
                nodalValues[i] = 0.0
                print p, nodalValues[i]
                cnt += 1
            else:
                nodalValues[i] = float(points[x])

        if cnt > 0:
            print '%i/%i of the grid points have been set to 0' % (cnt, gs.getSize())
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

#         err, _ = checkInterpolation(self.grid, alpha, nodalValues, epsilon=1e-8)
#
#         if len(err) > 0:
#             print "interpolation property not met"
#             pdb.set_trace()
        # -----------------------------------------

        return alpha

    def learnData(self, dataset=None, *args, **kws):
        # learning step
        self.notifyEventControllers(LearnerEvents.LEARNING_STARTED)

        if dataset is None:
            dataset = self.dataContainer

        # load data sets
        trainSubset = dataset.getTrainDataset()

        # learning step
        self.alpha = self.doLearningIteration(trainSubset)

        # update internal statistics
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
        self.trainErrors = self.evalError(trainSubset, alpha)

        # store interpolation quality
        if testSubset:
            self.testErrors = self.evalError(testSubset, alpha)

    def evalError(self, data, alpha):
        """
        computes the discrete L2-error of the sparse grid interpolant
        with respect to some MC samples at given time steps
        @param data: DataContainer samples
        @param alpha: numpy array hierarchical coefficients
        @return: mean squared error
        """
        # check if points are in [0, 1]^d
        allPoints = data.getPoints().array()
        allModelValues = data.getValues().array()

        points = np.ndarray((0, data.dim))
        modelValues = np.array([])
        for i, point in enumerate(allPoints):
            if np.all(point >= 0) and np.all(point <= 1):
                points = np.vstack((points, point))
                modelValues = np.append(modelValues, allModelValues[i])

        if allPoints.shape[1] == 0:
            return np.array([])

        surrogateValues = evalSGFunctionMulti(self.grid, alpha, points)

        # compute error
        return np.abs(surrogateValues - modelValues)

    def getSize(self, dtype="train"):
        if dtype == "train":
            return self.trainErrors.shape[0]
        else:
            return self.testErrors.shape[0]

    def getL2NormError(self, dtype="train"):
        """
        calculate L2-norm of error
        @return: last L2-norm of error
        """
        if dtype == "train":
            value = np.mean(self.trainErrors ** 2)
        else:
            value = np.mean(self.testErrors ** 2)
            
        return np.sqrt(value)

    def getL1NormError(self, dtype="train"):
        """
        calculate L2-norm of error
        @return: last L2-norm of error
        """
        if dtype == "train":
            value = np.mean(np.abs(self.trainErrors))
        else:
            value = np.mean(np.abs(self.testErrors))

        return value

    def getMaxError(self, dtype="train"):
        """
        calculate max error
        @return: max error
        """
        if dtype == "train":
            return self.trainErrors.max()
        else:
            return self.testErrors.max()

    def getMinError(self, dtype="train"):
        """
        calculate min error
        @return: min error
        """
        if dtype == "train":
            return self.trainErrors.min()
        else:
            return self.testErrors.min()

    def getErrorPerSample(self, dtype="train"):
        if dtype == "sample":
            return self.trainErrors
        else:
            return self.testErrors

    def getMSE(self, dtype="train"):
        if dtype == "train":
            return np.mean(self.trainErrors ** 2)
        else:
            return np.mean(self.testErrors ** 2)
