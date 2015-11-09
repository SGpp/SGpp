from pysgpp import (DataVector,
                    NaiveSampleGenerator,
                    LatinHypercubeSampleGenerator,
                    HaltonSampleGenerator,
                    SobolSampleGenerator,
                    ScrambledSobolSampleGenerator,
                    StratifiedSampleGenerator)
from Sample import Samples, SampleType, DistributionType
from Sampler import Sampler


class MCSampler(Sampler):

    def __init__(self, params, generator=NaiveSampleGenerator, *args, **kws):
        super(self.__class__, self).__init__()
        self.setParameters(params)
        self.__genCls = generator
        self.__genObj = generator(self._dim, *args, **kws)

    @classmethod
    def withNativeSampleGenerator(cls, params):
        return MCSampler(params, NaiveSampleGenerator)

#     @classmethod
#     def withStratifiedSampleGenerator(cls, params, strataPerDimension):
#         return MCSampler(params, StratifiedSampleGenerator, strataPerDimension)

    @classmethod
    def withLatinHypercubeSampleGenerator(cls, params, nSamples):
        return MCSampler(params, LatinHypercubeSampleGenerator, nSamples)

    @classmethod
    def withSobolSampleGenerator(cls, params):
        return MCSampler(params, SobolSampleGenerator)

    @classmethod
    def withHaltonSampleGenerator(cls, params):
        return MCSampler(params, HaltonSampleGenerator)

    @classmethod
    def withScrambledSobolSampleGenerator(cls, params):
        return MCSampler(params, ScrambledSobolSampleGenerator)

    def nextSamples(self, n=1):
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
