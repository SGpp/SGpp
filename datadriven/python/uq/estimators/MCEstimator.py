from Estimator import Estimator
import numpy as np
from scipy.stats import norm


class MCEstimator(Estimator):

    def __init__(self, npaths=100, percentile=1):
        """
        Constructor
        """
        Estimator.__init__(self)
        self.__npaths = npaths
        self.__percentile = percentile

    def getBootstrap(self, samples):
        ixs = np.random.randint(0, len(samples), len(samples))
        return samples[ixs]

    def mean(self, samples):
        """
        Compute the mean
        @param samples: numpy array
        """
        moments = np.ndarray(self.__npaths)
        for i in xrange(self.__npaths):
            bootstrap = self.getBootstrap(samples)
            moments[i] = np.mean(bootstrap)

        # error statistics
        if self.__npaths > 1:
            lower_percentile = np.percentile(moments, q=self.__percentile)
            upper_percentile = np.percentile(moments, q=100 - self.__percentile)
            err = max(lower_percentile, upper_percentile)
        else:
            err = lower_percentile = upper_percentile = np.Inf

        return {"value": np.mean(samples),
                "err": err,
                "confidence_interval": (lower_percentile, upper_percentile)}

    def var(self, samples):
        """
        Compute the variance
        @param samples: numpy array
        """
        moments = np.ndarray(self.__npaths)
        for i in xrange(self.__npaths):
            bootstrap = self.getBootstrap(samples)
            moments[i] = np.var(bootstrap)

        # error statistics
        if self.__npaths > 1:
            lower_percentile = np.percentile(moments, q=self.__percentile)
            upper_percentile = np.percentile(moments, q=100 - self.__percentile)
            err = max(lower_percentile, upper_percentile)
        else:
            err = lower_percentile = upper_percentile = np.Inf

        return {"value": np.var(samples),
                "err": err,
                "confidence_interval": (lower_percentile, upper_percentile)}
