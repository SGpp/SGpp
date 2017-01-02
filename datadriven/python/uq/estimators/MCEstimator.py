from Estimator import Estimator
import numpy as np
from scipy.stats import norm
from scikits.bootstrap import ci


class MCEstimator(Estimator):

    def __init__(self, npaths=100):
        """
        Constructor
        """
        Estimator.__init__(self)
        self._npaths = npaths

    def getBootstrap(self, samples):
        ixs = np.random.randint(0, len(samples), len(samples))
        return samples[ixs]

    def confidenceInterval(self, samples, alpha=0.01):
        if not np.all(np.abs(samples) < 1e-14):
            return ci(samples, alpha=alpha)
        else:
            return np.zeros(samples.shape)

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
