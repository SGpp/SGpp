from pysgpp_datadriven.uq.dists import (Dist, Uniform, Normal, TNormal, SGDEdist,
                          Lognormal, Beta, MultivariateNormal)
from pysgpp_datadriven.uq.transformation import (LinearTransformation,
                                   InverseCDFTransformation)

from DeterministicParameter import DeterministicParameter
from UncertainParameter import UncertainParameter
from pysgpp_datadriven.uq.transformation.JointTransformation import JointTransformation


class ParameterDescriptor(object):

    def __init__(self):
        self._name = None
        self._value = None
        self._correlatedTo = None

    def isCalled(self, name):
        if ',' in name:
            self._name = [x.strip() for x in name.split(',')]
        else:
            self._name = name
        return self

    def hasValue(self, value):
        self._value = value
        return self

    def isCorrelatedTo(self, a):
        self._correlatedTo = a
        return self

    def getCorrelations(self):
        return self._correlatedTo

    def andGetResult(self):
        raise NotImplementedError()


class DeterministicParameterDescriptor(ParameterDescriptor):
    """
    Descriptor, specifying one deterministic parameter
    """

    def __init__(self):
        super(DeterministicParameterDescriptor, self).__init__()

    def andGetResult(self):
        trans = LinearTransformation(self._value, self._value)
        return DeterministicParameter(self._name, self._value, trans)


class UncertainParameterDesciptor(ParameterDescriptor):
    """
    Descriptor, specifying one uncertain parameter
    """

    def __init__(self):
        super(UncertainParameterDesciptor, self).__init__()
        self._dist = None
        self.__trans = None

    def withSGDEDistribution(self, dist):
        if isinstance(dist, SGDEdist):
            self._dist = dist
        else:
            raise AttributeError('the argument needs to be a SGDEDist object')
        return self

    def withSGDEConfig(self, config, *args, **kws):
        """
        Estimates the density from training data
        @param config: configuration file for density estimation
        @return: self
        """
        self._dist = SGDEdist(config, *args, **kws)
        return self

    def withUniformDistribution(self, a, b):
        self._dist = Uniform(a, b)
        return self

    def withTNormalDistribution(self, mu, sigma, a, b):
        self._dist = TNormal(mu, sigma, a, b)
        return self

    def withNormalDistribution(self, mu, sigma, alpha):
        self._dist = Normal.by_alpha(mu, sigma, alpha)
        return self

    def withLognormalDistribution(self, mu, sigma, alpha):
        self._dist = Lognormal.by_alpha(mu, sigma, alpha)
        return self

    def withBetaDistribution(self, p, q, accLevel=0., width=1.):
        self._dist = Beta(p, q, accLevel, width)
        return self

    def withMultivariateNormalDistribution(self, mu, cov, a, b):
        self._dist = MultivariateNormal(mu, cov, a, b)
        return self

    def withDistribution(self, dist):
        if issubclass(dist.__class__, Dist):
            self._dist = dist
            return self
        else:
            raise TypeError('dist has to have type dists.Dist')

    def withLinearTransformation(self):
        if self._dist is not None:
            a, b = self._dist.getBounds()
            self.__trans = LinearTransformation(a, b)
        else:
            raise AttributeError('the distribution of "%s" is not specified \
                                  yet but it is needed to know to apply the \
                                  LinearTransformation' % self.__name)
        return self

    def withInverseCDFTransformation(self):
        if self._dist is not None:
            self.__trans = InverseCDFTransformation(self._dist)
        else:
            raise AttributeError('the distribution of "%s" is not specified \
                                  yet but it is needed to know to apply the \
                                  InverseCDFTransformation' % self.__name)
        return self

    def andGetResult(self):
        if self._dist is None:
            raise Exception("No distribution specified for parameter '%s'" % self._name)
        # if there is no transformation defined, use the
        # linear transformation as standard
        if self.__trans is None:
            bounds = self._dist.getBounds()
            if self._dist.getDim() == 1:
                self.__trans = LinearTransformation(bounds[0], bounds[1])
            else:
                self.__trans = JointTransformation()
                for i in xrange(self._dist.getDim()):
                    a, b = bounds[i]
                    self.__trans.add(LinearTransformation(a, b))

        # check if there are enough identifiers given
        if self._dist.getDim() > 1 and self._dist.getDim() != len(self._name):
            raise AttributeError("not enough names provided; given %i (%s), expected %i" % (len(self._name), ", ".join(self._name), self._dist.getDim()))

        # check if there is and SGDE involved that need this
        # transformation
        if isinstance(self._dist, SGDEdist):
            self._dist.transformation = self.__trans

        return UncertainParameter(self._name, self._dist,
                                  self.__trans, self._value)
