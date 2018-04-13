from pysgpp.extensions.datadriven.uq.transformation import RosenblattTransformation, \
    JointTransformation
from pysgpp.extensions.datadriven.uq.dists import Uniform, J
from pysgpp.extensions.datadriven.uq.transformation.LinearTransformation import LinearTransformation


class SparseGridEstimationStrategy(object):

    def _extractPDFforMomentEstimation(self, U, T):
        dists = []
        jointTrans = []
        vol = 1.
        # check if importance sampling has been used for some parameters
        for i, trans in enumerate(T.getTransformations()):
            # if this is the case replace them by a uniform distribution
            if isinstance(trans, RosenblattTransformation):
                for _ in xrange(trans.getSize()):
                    dists.append(Uniform(0, 1))
                    jointTrans.append(LinearTransformation(0.0, 1.0))
            else:
                vol *= trans.vol()
                dists.append(U.getDistributions()[i])
                jointTrans.append(trans)
        return vol, J(dists), jointTrans

    def mean(self, grid, alpha, U, T, *args, **kws):
        """
        Compute the expectation value

        @param grid: Grid
        @param alpha: DataVector coefficients
        @param U: Dist distribution
        @param T: Transformation function
        """
        raise NotImplementedError()

    def var(self, grid, alpha, U, T, mean, *args, **kws):
        """
        Compute the variance

        @param grid: Grid
        @param alpha: DataVector coefficients
        @param U: Dist distribution
        @param T: Transformation function
        @param mean: double mean of sparse grid function
        """
        raise NotImplementedError()

    def toJson(self):
        """
        Returns a string that represents the object
        """
        return '{"module" : "' + self.__module__ + '"}'

    @classmethod
    def fromJson(cls, jsonObject):
        import pysgpp.extensions.datadriven.uq.estimators as estimators
        if jsonObject['module'] == 'estimators.IntegralStrategy':
            return estimators.IntegralStrategy()
        elif jsonObject['module'] == 'estimators.MonteCarloStrategy':
            return estimators.MonteCarloStrategy()
        elif jsonObject['module'] == 'estimators.CollocationPointsStrategy':
            return estimators.CollocationPointsStrategy()
        else:
            raise TypeError('Unknown estimation strategy "%s" => Please \
                            register it in fromJson function' %
                            jsonObject['module'])
