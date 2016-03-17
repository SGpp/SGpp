import numpy as np

from pysgpp.extensions.datadriven.uq.transformation.LinearTransformation import LinearTransformation
from pysgpp.extensions.datadriven.uq.transformation.JointTransformation import JointTransformation
from pysgpp.extensions.datadriven.uq.dists.Dist import Dist
from pysgpp.extensions.datadriven.uq.operations.general import isList

class EstimatedDist(Dist):

    def __init__(self, trainData):
        super(EstimatedDist, self).__init__()

        if isList(trainData) or len(trainData.shape) == 1:
            trainData = np.array([trainData]).reshape(len(trainData), 1)

        self.trainData = trainData
        self.dim = trainData.shape[1]
        self.trans = self.computeLinearTransformation(trainData)

    def rvs(self, n=1):
        unif = np.random.rand(self.dim * n).reshape(n, self.dim)
        return self.ppf(DataMatrix(unif))

    @classmethod
    def computeLinearTransformation(self, trainData):
        num_dims = trainData.shape[1]
        bounds = np.ndarray((num_dims, 2))
        bounds[:, 0] = trainData.min(axis=0) * -0.95
        bounds[:, 1] = trainData.max(axis=0) * 1.05

        # init linear transformation
        trans = JointTransformation()
        for idim in xrange(num_dims):
            trans.add(LinearTransformation(bounds[idim, 0], bounds[idim, 1]))

        return trans

    def getDistributions(self):
        return [self]

    def getBounds(self):
        return self.trans.getBounds()

    def getDim(self):
        return self.dim

    def getSamples(self):
        return self.trainData

    def __str__(self):
        return "EstimatedDist"
