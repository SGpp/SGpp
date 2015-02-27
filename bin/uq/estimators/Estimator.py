from MonteCarloStrategy import MonteCarloStrategy


class Estimator(object):

    def __init__(self, params, strategy=None):
        """
        Constructor
        @param params: parameter set
        @param k: moments to be estimated {E(f^k) | k \in [1, ..., k]}
        @param strategy: estimation strategy
        """
        self.__params = params.activeParams()
        self.__T = self.__params.getJointTransformation()
        self.__vol = self.__T.vol()
        self.__U = self.__params.getIndependentJointDistribution()
        if strategy is None:
            self.__strategy = MonteCarloStrategy()
        else:
            self.__strategy = strategy

    def mean(self, grid, alpha, *args, **kws):
        """
        Estimate moments of the given sparse grid function
        @param grid: Grid
        @param alpha: coefficients vector
        @return: tuple(moment, error)
        """
        return self.__strategy.mean(self.__vol, grid, alpha,
                                    self.__U, self.__T,
                                    *args, **kws)

    def var(self, grid, alpha, *args, **kws):
        """
        Estimate moments of the given sparse grid function
        @param grid: Grid
        @param alpha: coefficients vector
        @return: tuple(moment, error)
        """
        return self.__strategy.mean(self.__vol, grid, alpha,
                                    self.__U, self.__T,
                                    *args, **kws)

    def setStrategy(self, strategy):
        self.__strategy = strategy

    def getStrategy(self):
        return self.__strategy

    def setNumberOfMoments(self, k):
        self.__k = k

    def getNumberOfMoments(self):
        return self.__k
