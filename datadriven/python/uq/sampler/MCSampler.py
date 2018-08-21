from pysgpp import (DataVector,
                    NaiveSampleGenerator,
                    LatinHypercubeSampleGenerator,
                    HaltonSampleGenerator,
                    StratifiedSampleGenerator)
from Sample import Samples, SampleType, DistributionType
from Sampler import Sampler


class MCSampler(Sampler):

    def __init__(self, params=None, samples=None, generator=NaiveSampleGenerator, *args, **kws):
        super(self.__class__, self).__init__()
        self.samples = samples
        self.setParameters(params)

        self.__genCls = generator
        if generator != "numpy":
            self.__genObj = generator(self._dim, *args, **kws)

    @classmethod
    def withNaiveSampleGenerator(cls, params):
        return MCSampler(params, None, NaiveSampleGenerator)

    @classmethod
    def withStratifiedSampleGenerator(cls, params, strataPerDimension):
        return MCSampler(params, None, StratifiedSampleGenerator, strataPerDimension)

    @classmethod
    def withLatinHypercubeSampleGenerator(cls, params, nSamples):
        return MCSampler(params, None, LatinHypercubeSampleGenerator, nSamples)

    @classmethod
    def withHaltonSampleGenerator(cls, params):
        return MCSampler(params, None, HaltonSampleGenerator)

    @classmethod
    def withNumpySampleGenerator(cls, params):
        return MCSampler(params, None, generator="numpy")

    def nextSamples(self, n=1):
        if self.__genCls == "numpy":
            ans = Samples(self._params, dtype=DistributionType.PROBABILISTICDIST)
            dist = self._params.activeParams().getIndependentJointDistribution()
            samples = dist.rvs(n)
            for sample in samples:
                ans.add(sample, dtype=SampleType.ACTIVEPROBABILISTIC)
        else:
            p = DataVector(self._dim)
            ans = Samples(self._params, dtype=DistributionType.UNITUNIFORM)
            U = self._params.activeParams().getIndependentJointDistribution()
            for _ in xrange(n):
                self.__genObj.getSample(p)

                # transform it to the probabilistic space
                q = U.ppf(p.array())
                # add it to the output
                ans.add(q, dtype=SampleType.ACTIVEPROBABILISTIC)

        return ans

    def hasMoreSamples(self):
        return True

    def learnData(self, *args, **kws):
        return
