// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/NormStrategy.hpp>
#include <sgpp/combigrid/algebraic/FirstMomentNormStrategy.hpp>
#include <sgpp/combigrid/algebraic/SecondMomentNormStrategy.hpp>
#include <sgpp/combigrid/algebraic/VarianceNormStrategy.hpp>
#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class VarianceNormStrategy : public NormStrategy<FloatTensorVector> {
 public:
  VarianceNormStrategy(std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction,
                       size_t numDims, sgpp::combigrid::SingleFunction weightFunction,
                       bool isOrthogonal,
                       sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  VarianceNormStrategy(std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction,
                       sgpp::combigrid::WeightFunctionsCollection& weightFunctions,
                       bool isOrthogonal,
                       sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  VarianceNormStrategy(sgpp::combigrid::OrthogonalBasisFunctionsCollection& basisFunctions,
                       sgpp::combigrid::WeightFunctionsCollection& weightFunctions,
                       bool isOrthogonal,
                       sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  ~VarianceNormStrategy() override;

  /**
   * computes the variance of the given tensor vector
   *
   * @param vector tensor
   * @return variance
   */
  double norm(FloatTensorVector& vector) override;

 private:
  FirstMomentNormStrategy firstMoment;
  SecondMomentNormStrategy secondMoment;
};

} /* namespace combigrid */
} /* namespace sgpp */
