import numpy as np
import pysgpp.extensions.datadriven.uq.jsonLib as ju
from copy import copy


class SampleType:
    """
    Describes the type of a sample:
    An active sample is a sample that just contains uncertain
    parameters for which no deterministic value has been specified.
    This allows to increase the number of uncertain parameters step
    by step and reuse the results. Because the sparse grid is used
    just on the unit hyper cube we distinguish four different
    classes of samples:
    1) active unit: the tuple just contains the active parameters,
                    which are scaled to the unit hyper cube
    2) active probabilistic: the tuple contains the active parameters,
                             which are scaled to the probabilistic space
    3) expanded unit: the tuple contains all parameters, and they are
                      scaled to the unit hyper cube
    4) expanded probabilistic: the tuple contains all parameters, and they
                               are scaled to the probabilistic space
    """
    ACTIVEUNIT = 1
    ACTIVEPROBABILISTIC = 2
    EXPANDEDUNIT = 3
    EXPANDEDPROBABILISTIC = 4


class DistributionType:
    """
    Describes how the samples are drawn:
    1) unit and uniform: the samples are all in the unit hyper cube and are
                         sampled uniformly
    2) probabilistic and uniform: the samples are in the probabilistic space
                                  but sampled uniformly
    3) probabilistic and distributed: the samples are in the probabilistic
                                      space and sampled according to their
                                      distribution

    Examples:
    - MCSampler draws samples in (1)
    - rvs() draws samples in (3)
    - after applying transformation to samples in (3) they can either be in
        (2) if the transformation is linear
        (1) if the transformation is the inverse CDF transformation
    """
    UNITUNIFORM = 1
    PROBABILISTICUNIFORM = 2
    PROBABILISTICDIST = 3


class Samples(object):

    def __init__(self, params, dtype=DistributionType.UNITUNIFORM):
        self._samples = {}

        self._params = params
        self._dtype = dtype
        self._isUnit = True
        self._isActive = True
        self._dim = None
        self._dim = params.getStochasticDim()

    def combine(self, samples):
        for sample in samples:
            self.addSample(sample)

    def add(self, p, dtype=SampleType.ACTIVEUNIT):
        sample = Sample(self._params, tuple(p), dtype=dtype)
        self.addSample(sample)

    def removeSet(self, samples):
        for sample in samples:
            self.removeSample(sample)

    def addSample(self, sample):
        self._samples[sample] = sample

    def removeSample(self, sample):
        if sample in self._samples:
            del self._samples[sample]

    def selectActiveElements(self):
        self._isActive = True
        return self

    def selectUncertainElements(self):
        self._isActive = False
        return self

    def selectUnitSpace(self):
        self._isUnit = True
        return self

    def selectProbabilisticSpace(self):
        self._isUnit = False
        return self

    def __getitem__(self, i):
        keys, values = zip(*self._samples.items())
        if isinstance(i, slice):
            # copy object but samples
            ans = Samples(self._params, dtype=self._dtype)
            ans._isUnit = self._isUnit
            ans._isActive = self._isActive
            # copy samples
            ans._samples = {}
            for key, value in zip(keys[i], values[i]):
                ans._samples[key] = value
            return ans
        else:
            return values[i]

    def __len__(self):
        return len(self._samples)

    def getDim(self, *args, **kws):
        self._samples.itervalues().next().getDim(*args, **kws)

    def array(self):
        self._samples

    def ndarray(self):
        samples = np.ndarray([len(self._samples), self._dim], dtype='float')
        if self._isUnit:
            if self._isActive:
                # return just active tuples from the unit hypercube
                for i, sample in enumerate(self._samples.values()):
                    samples[i, :] = sample.getActiveUnit()
            else:
                # return expanded tuples from the unit hypercube
                for i, sample in enumerate(self._samples.values()):
                    samples[i, :] = sample.getExpandedUnit()
        else:
            if self._isActive:
                # return just active tuples from the probabilitic space
                for i, sample in enumerate(self._samples.values()):
                    samples[i, :] = sample.getActiveProbabilistic()
            else:
                # return expanded tuples from the probabilistic space
                for i, sample in enumerate(self._samples.values()):
                    samples[i, :] = sample.getExpandedProbabilistic()

        return samples

    def __iter__(self):
        return SamplesIterator(self._samples)

    def __str__(self):
        return "SampleSet: nsamples = %i, ndims = %i" % (len(self._samples),
                                                         self._dim)

class SamplesIterator(object):
    """
    Iterator class
    """
    def __init__(self, samples):
        self.samples = samples
        self.__current = 0
        self.__values = self.samples.values()

    def next(self):
        if self.__current == len(self.samples):
            raise StopIteration()

        ans = self.__values[self.__current]
        self.__current += 1
        return ans


class Sample(object):

    def __init__(self, params=None, sample=None,
                 dtype=SampleType.ACTIVEPROBABILISTIC):
        """
        constructor
        @param params: ParameterSet
        @param sample: numpy array, tuple or list
        @param dtype: SampleType
        """
        self.__activeUnit = None
        self.__activeProb = None
        self.__expandedUnit = None
        self.__expandedProb = None
        # load the tuples if they can be
        if params is not None and sample is not None:
            trans = params.activeParams().getJointTransformation()
            if dtype == SampleType.ACTIVEUNIT:
                # sample contains just active dimensions in [0, 1]^d
                self.__activeUnit = np.array(sample, dtype='float')
                self.__activeProb = trans.unitToProbabilistic(sample)
                self.__expandedUnit = params.expandUnitParameter(self.__activeUnit)
                self.__expandedProb = params.expandProbabilisticParameter(self.__activeProb)
            elif dtype == SampleType.ACTIVEPROBABILISTIC:
                # sample contains just active dimensions in \Gamma
                self.__activeProb = np.array(sample, dtype='float')
                self.__activeUnit = trans.probabilisticToUnit(sample)
                self.__expandedUnit = params.expandUnitParameter(self.__activeUnit)
                self.__expandedProb = params.expandProbabilisticParameter(self.__activeProb)
            elif dtype == SampleType.EXPANDEDUNIT:
                # sample contains all dimensions in [0, 1]^d
                self.__expandedUnit = np.array(sample, dtype='float')
                self.__activeUnit = params.extractActiveTuple(sample)
                self.__activeProb = trans.unitToProbabilistic(self.__activeUnit)
                self.__expandedProb = params.expandProbabilisticParameter(self.__activeProb)
            elif dtype == SampleType.EXPANDEDPROBABILISTIC:
                # sample contains all dimensions in \Gamma
                self.__expandedProb = np.array(sample, dtype='float')
                self.__activeProb = params.extractActiveTuple(sample)
                self.__activeUnit = trans.probabilisticToUnit(self.__activeProb)
                self.__expandedUnit = params.expandUnitParameter(self.__activeUnit)
            else:
                raise AttributeError('dtype "%s" is not known for a \
                                      Sample' % dtype)

    def init(self, activeUnit, activeProb, expandedUnit, expandedProb):
        self.__activeProb = activeProb
        self.__activeUnit = activeUnit
        self.__expandedProb = expandedProb
        self.__expandedUnit = expandedUnit

    def getActiveUnit(self):
        return self.__activeUnit

    def getActiveProbabilistic(self):
        return self.__activeProb

    def getExpandedUnit(self):
        return self.__expandedUnit

    def getExpandedProbabilistic(self):
        return self.__expandedProb

    def getStochasticDim(self):
        return len(self.__activeUnit)

    def getExpandedDim(self):
        return len(self.__expandedUnit)

    def getDim(self, dtype=SampleType.ACTIVEUNIT):
        if dtype in [SampleType.ACTIVEPROBABILISTIC,
                     SampleType.ACTIVEUNIT]:
            return self.getStochasticDim()
        elif dtype in [SampleType.EXPANDEDPROBABILISTIC,
                       SampleType.EXPANDEDUNIT]:
            return self.getExpandedDim()
        else:
            raise AttributeError('the given SampleType %i is not known' %
                                 dtype)

    def getValue(self, dtype=SampleType.ACTIVEUNIT):
        if dtype == SampleType.ACTIVEUNIT:
            return self.getActiveUnit()
        elif dtype == SampleType.ACTIVEPROBABILISTIC:
            return self.getActiveProbabilistic()
        elif dtype == SampleType.EXPANDEDUNIT:
            return self.getExpandedUnit()
        elif dtype == SampleType.EXPANDEDPROBABILISTIC:
            return self.getExpandedProbabilistic()
        else:
            raise AttributeError('the given SampleType %i is not known' %
                                 dtype)

    def __str__(self):
        return str(self.getActiveUnit())

    def __hash__(self):
        return tuple(self.getExpandedUnit()).__hash__()

    def __eq__(self, sample):
        return hash(self) == hash(sample)

    # ---------------------------------------------------------------------------
    # Json serialization methods
    # ---------------------------------------------------------------------------
    def toJson(self):
        """
        Returns a string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        # serialize
        for attrName in ('_Sample__activeUnit',
                         '_Sample__activeProb',
                         '_Sample__expandedUnit',
                         '_Sample__expandedProb'):
            attrValue = self.__getattribute__(attrName)
            serializationString += ju.parseAttribute(attrValue, attrName)

        return "{" + serializationString.rstrip(",\n") + "}"

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the Sample object from the json object
        with its attributes.
        @param jsonObject: json object
        @return: the restored Sample object
        """
        sample = Sample()
        key = '_Sample__activeUnit'
        if key in jsonObject:
            activeUnit = np.array(jsonObject[key], dtype='float')
        key = '_Sample__activeProb'
        if key in jsonObject:
            activeProb = np.array(jsonObject[key], dtype='float')
        key = '_Sample__expandedUnit'
        if key in jsonObject:
            expandedUnit = np.array(jsonObject[key], dtype='float')
        key = '_Sample__expandedProb'
        if key in jsonObject:
            expandedProb = np.array(jsonObject[key], dtype='float')

        # initialize the sample
        sample.init(activeUnit, activeProb, expandedUnit, expandedProb)
        return sample
