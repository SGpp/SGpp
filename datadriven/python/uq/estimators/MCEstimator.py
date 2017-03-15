from Estimator import Estimator
import numpy as np
from scipy.stats import norm


class MCEstimator(Estimator):

    def __init__(self, npaths=10):
        """
        Constructor
        @param n: number of samples
        """
        Estimator.__init__(self)
        self._npaths = npaths
#         self._c = norm.ppf(1. - beta / 2.)

    def getBootstrap(self, samples):
        ixs = np.random.randint(0, len(samples), len(samples))
        return samples[ixs]

    def mean(self, samples):
        """
        Compute the mean
        @param samples: numpy array
        """
        moments = np.ndarray(self._npaths)
        for i in xrange(self._npaths):
            moments[i] = np.mean(self.getBootstrap(samples))

        # error statistics
        if self._npaths > 1:
            err = np.var(moments, ddof=1)
        else:
            err = np.Inf

        return np.mean(moments), err

    def var(self, samples):
        moments = np.ndarray(self._npaths)
        for i in xrange(self._npaths):
            moments[i] = np.var(self.getBootstrap(samples), ddof=1)

        # error statistics
        if self._npaths > 1:
            err = np.var(moments, ddof=1)
        else:
            err = np.Inf

        return np.mean(moments), err
