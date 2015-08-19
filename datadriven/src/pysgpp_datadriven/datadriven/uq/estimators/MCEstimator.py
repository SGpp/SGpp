from Estimator import Estimator
import numpy as np
from scipy.stats import norm


class MCEstimator(Estimator):

    def __init__(self,
                 n=1000,
                 npaths=10,
                 isPositive=False,
                 beta=0.01):
        """
        Constructor
        @param params: parameter set
        """
        Estimator.__init__(self)
        self._npaths = npaths
        self._n = n
        self._isPositive = isPositive

        self._c = norm.ppf(1. - beta / 2.)

    def getSubset(self, samples):
        np.random.shuffle(samples)
        return samples[:self._n]

    def mean(self, samples):
        """
        Compute the mean
        @param samples: numpy array
        """
        moments = np.ndarray(self._npaths)
        for i in xrange(self._npaths):
            moments[i] = np.mean(self.getSubset(samples))

        # error statistics
        if self._npaths > 1:
            err_clt = self._c * np.sqrt(np.var(moments, ddof=1) / len(moments))
        else:
            err_clt = np.Inf

        return np.mean(moments), err_clt

    def var(self, samples):
        moments = np.ndarray(self._npaths)
        for i in xrange(self._npaths):
            moments[i] = np.var(self.getSubset(samples), ddof=1)

        # error statistics -> not available
        err_clt = np.Inf

        return np.mean(moments), err_clt
