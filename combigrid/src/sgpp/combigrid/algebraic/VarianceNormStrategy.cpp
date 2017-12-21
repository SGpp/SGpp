// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/algebraic/VarianceNormStrategy.hpp>
#include <sgpp/combigrid/algebraic/FirstMomentNormStrategy.hpp>
#include <sgpp/combigrid/algebraic/SecondMomentNormStrategy.hpp>
#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

VarianceNormStrategy::VarianceNormStrategy(
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction, size_t numDims,
    sgpp::combigrid::SingleFunction weightFunction, bool isOrthogonal,
    sgpp::base::DataVector const& bounds)
    : firstMoment(basisFunction, numDims, weightFunction, isOrthogonal, bounds),
      secondMoment(basisFunction, numDims, weightFunction, isOrthogonal, bounds) {}

VarianceNormStrategy::VarianceNormStrategy(
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction,
    sgpp::combigrid::WeightFunctionsCollection& weightFunctions, bool isOrthogonal,
    sgpp::base::DataVector const& bounds)
    : firstMoment(basisFunction, weightFunctions, isOrthogonal, bounds),
      secondMoment(basisFunction, weightFunctions, isOrthogonal, bounds) {}

VarianceNormStrategy::VarianceNormStrategy(
    sgpp::combigrid::OrthogonalBasisFunctionsCollection& basisFunctions,
    sgpp::combigrid::WeightFunctionsCollection& weightFunctions, bool isOrthogonal,
    sgpp::base::DataVector const& bounds)
    : firstMoment(basisFunctions, weightFunctions, isOrthogonal, bounds),
      secondMoment(basisFunctions, weightFunctions, isOrthogonal, bounds) {}

VarianceNormStrategy::~VarianceNormStrategy() {}

double VarianceNormStrategy::norm(FloatTensorVector& vector) {
  double second = secondMoment.norm(vector);
  double first = firstMoment.norm(vector);

  return second - first * first;
}

} /* namespace combigrid */
} /* namespace sgpp */
