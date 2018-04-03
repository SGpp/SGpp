import numpy as np

from pysgpp.extensions.datadriven.uq.transformation.LinearTransformation import LinearTransformation
from pysgpp.extensions.datadriven.uq.transformation.JointTransformation import JointTransformation
from pysgpp.extensions.datadriven.uq.dists.Dist import Dist
from pysgpp.extensions.datadriven.uq.operations.general import isList, \
    isNumerical

class EstimatedDist(Dist):

    def __init__(self, numDims, trainData=None, bounds=None):
        super(EstimatedDist, self).__init__()

        self.trainData = None
        self.dim = numDims
        self.trans = None
        
        if trainData is not None:
            if isList(trainData) or len(trainData.shape) == 1:
                trainData = np.array([trainData]).reshape(len(trainData), 1)

            self.trainData = trainData
            self.dim = trainData.shape[1]

            if bounds is None:
                # estimate bounds from data
                bounds = np.vstack((np.min(trainData, axis=0), np.max(trainData, axis=0))).T

        if bounds is None:
            self.bounds = np.array([[0, 1]] * self.dim, dtype="float")
        else:
            self.bounds = bounds

        self.trans = self.computeLinearTransformation(self.bounds)

    def rvs(self, n=1):
        unif = np.random.rand(self.dim * n).reshape(n, self.dim)
        return self.ppf(unif)

    @classmethod
    def computeLinearTransformation(self, bounds):
        if len(bounds.shape) == 1:
            bounds = np.array([bounds])
        if bounds.shape[1] != 2:
            raise AttributeError("EstimatedDist - bounds have the wrong shape")

        # init linear transformation
        trans = JointTransformation()
        for idim in xrange(bounds.shape[0]):
            trans.add(LinearTransformation(bounds[idim, 0], bounds[idim, 1]))

        return trans

    def _convertEvalPoint(self, x):
        # convert the parameter to the right format
        if self.dim == 1:
            if isNumerical(x):
                x = np.array([[x]])
            elif isList(x) or len(x.shape) == 1:
                x = np.array([x]).reshape(len(x), 1)
        else:
            x = np.array(x)
            if len(x.shape) == 1:
                x = np.array([x])
        return x

    def marginalizeToDimX(self):
        raise NotImplementedError()

    def getDistributions(self):
        return [self]

    def getBounds(self):
        if self.getDim() == 1 and len(self.bounds.shape) > 1:
            return self.bounds[0]
        else:
            return self.bounds

    def getDim(self):
        return self.dim

    def getSamples(self):
        return self.trainData

    def __str__(self):
        return "EstimatedDist"
