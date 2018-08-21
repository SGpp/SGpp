import numpy as np

from pysgpp.extensions.datadriven.uq.dists import (Dist, Uniform, Normal, TNormal, SGDEdist,
                          Lognormal, Beta, MultivariateNormal, TLognormal)
from pysgpp.extensions.datadriven.uq.transformation import (LinearTransformation,
                                                            RosenblattTransformation)

from DeterministicParameter import DeterministicParameter
from UncertainParameter import UncertainParameter
from pysgpp.extensions.datadriven.uq.transformation.JointTransformation import JointTransformation
from pysgpp.extensions.datadriven.uq.dists.DataDist import DataDist
from pysgpp.extensions.datadriven.uq.sampler.Sample import DistributionType
from pysgpp.pysgpp_swig import OrthogonalPolynomialBasis1DConfiguration, \
    OrthogonalPolynomialBasisType_HERMITE, \
    OrthogonalPolynomialBasisType_JACOBI, \
    OrthogonalPolynomialBasisType_LEGENDRE, \
    OrthogonalPolynomialBasisType_BOUNDED_LOGNORMAL, \
    OrthogonalPolynomialBasis1D


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

    def withUniformDistribution(self, a, b):
        self._dist = Uniform(a, b)
        return self

    def withTNormalDistribution(self, mu, sigma, a, b):
        self._dist = TNormal(mu, sigma, a, b)
        return self

    def withNormalDistribution(self, mu, sigma, alpha):
        self._dist = Normal.by_alpha(mu, sigma, alpha)
        return self

    def withTLognormalDistribution(self, mu, sigma, alpha):
        self._dist = TLognormal.by_alpha(mu, sigma, alpha)
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

    def withSampleDistribution(self, samples):
        self._dist = DataDist(samples)
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

    def withRosenblattTransformation(self):
        if self._dist is not None:
            self.__trans = RosenblattTransformation(self._dist)
        else:
            raise AttributeError('the distribution of "%s" is not specified yet but it is needed to know to apply the Rosenblatt transformation' % self.__name)
        return self

    def withTransformation(self, trans):
        self.__trans = trans
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

        # check wiener askey scheme
        partOfWienerAskeyScheme = False
        orthogPolyConfig = OrthogonalPolynomialBasis1DConfiguration()
        if isinstance(self.__trans, RosenblattTransformation) or \
                isinstance(self._dist, Uniform):
            orthogPolyConfig.polyParameters.type_ = OrthogonalPolynomialBasisType_LEGENDRE
            orthogPolyConfig.polyParameters.lowerBound_ = self._dist.getBounds()[0]
            orthogPolyConfig.polyParameters.upperBound_ = self._dist.getBounds()[1]
            orthogPoly = OrthogonalPolynomialBasis1D(orthogPolyConfig)
        elif isinstance(self._dist, Normal):
            orthogPolyConfig.polyParameters.type_ = OrthogonalPolynomialBasisType_HERMITE
            orthogPoly = OrthogonalPolynomialBasis1D(orthogPolyConfig)
        elif isinstance(self._dist, Beta):
            orthogPolyConfig.polyParameters.type_ = OrthogonalPolynomialBasisType_JACOBI
            orthogPolyConfig.polyParameters.alpha_ = self._dist.alpha()
            orthogPolyConfig.polyParameters.beta_ = self._dist.beta()
            orthogPoly = OrthogonalPolynomialBasis1D(orthogPolyConfig)
        elif isinstance(self._dist, TLognormal):
            orthogPolyConfig.polyParameters.type_ = OrthogonalPolynomialBasisType_BOUNDED_LOGNORMAL
            orthogPolyConfig.polyParameters.logmean_ = np.log(self._dist.mu)
            orthogPolyConfig.polyParameters.stddev_ = self._dist.sigma
            orthogPolyConfig.polyParameters.lowerBound_ = self._dist.getBounds()[0]
            orthogPolyConfig.polyParameters.upperBound_ = self._dist.getBounds()[1]
            orthogPoly = OrthogonalPolynomialBasis1D(orthogPolyConfig)
        else:
            orthogPoly = None

        return UncertainParameter(self._name, self._dist,
                                  self.__trans, self._value,
                                  orthogPoly)
