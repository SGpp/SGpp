
class Estimator(object):

    def mean(self, *args, **kws):
        """
        Estimate the mean
        @return: tuple(moment, error)
        """
        raise NotImplementedError()

    def var(self, *args, **kws):
        """
        Estimate the variance
        """
        raise NotImplementedError()
