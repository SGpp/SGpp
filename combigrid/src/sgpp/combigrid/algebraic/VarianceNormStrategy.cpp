// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/algebraic/VarianceNormStrategy.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

VarianceNormStrategy::VarianceNormStrategy(
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction, size_t numDims,
    bool isOrthogonal)
    : firstMoment(basisFunction, numDims, isOrthogonal),
      secondMoment(basisFunction, numDims, isOrthogonal) {}
VarianceNormStrategy::VarianceNormStrategy(
    std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& basisFunctions,
    bool isOrthogonal)
    : firstMoment(basisFunctions, isOrthogonal), secondMoment(basisFunctions, isOrthogonal) {}

VarianceNormStrategy::~VarianceNormStrategy() {}

} /* namespace combigrid */
} /* namespace sgpp */
