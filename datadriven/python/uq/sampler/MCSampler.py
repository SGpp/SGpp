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

    def __init__(self, params=None, samples=None, generator=NaiveSampleGenerator, *args, **kws):
        super(self.__class__, self).__init__()
        self.samples = samples

        if samples is None:
            self.setParameters(params)
            self.__genCls = generator
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
    def withSobolSampleGenerator(cls, params):
        return MCSampler(params, None, SobolSampleGenerator)

    @classmethod
    def withHaltonSampleGenerator(cls, params):
        return MCSampler(params, None, HaltonSampleGenerator)

    @classmethod
    def withScrambledSobolSampleGenerator(cls, params):
        return MCSampler(params, None, ScrambledSobolSampleGenerator)

    @classmethod
    def withSamples(cls, samples):
        return MCSampler(samples=samples)

    def nextSamples(self, n=1):
        if self.samples is not None:
            # use bootstrapping to generate data
#             ixs = np.random.randint(0, self.samples.shape[0], n)
            ans = Samples(None, dtype=DistributionType.PROBABILISTICUNIFORM)
            for i in xrange(self.samples.shape[0]):
                p = self.samples[i, :]
                ans.add(p, dtype=SampleType.ACTIVEUNIT)

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
