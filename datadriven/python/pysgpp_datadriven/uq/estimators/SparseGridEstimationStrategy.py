from pysgpp.extensions.datadriven.uq.transformation import InverseCDFTransformation
from pysgpp.extensions.datadriven.uq.dists import Uniform, J


class SparseGridEstimationStrategy(object):

    def _extractPDFforMomentEstimation(self, U, T):
        dists = U.getDistributions()
        vol = 1.
        # check if importance sampling has been used for some parameters
        for i, trans in enumerate(T.getTransformations()):
            # if this is the case replace them by a uniform distribution
            if isinstance(trans, InverseCDFTransformation):
                dists[i] = Uniform(0, 1)
            else:
                vol *= trans.vol()
        return vol, J(dists)

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
