// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/NormStrategy.hpp>
#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class FirstMomentNormStrategy : public NormStrategy<FloatTensorVector> {
 public:
  FirstMomentNormStrategy(
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction, size_t numDims,
      bool isOrthogonal);
  FirstMomentNormStrategy(
      std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& basisFunctions,
      bool isOrthogonal);

  virtual ~FirstMomentNormStrategy();

  double norm(FloatTensorVector& vector) {
    if (isOrthogonal) {
      auto values = vector.getValues();
      MultiIndex ix(values->getNumDimensions(), 0);
      return values->get(ix).value();
    } else {
      return computeMean(vector);
    }
  }

 private:
  bool isOrthogonal;
  std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>> basisFunctions;

  double quad(MultiIndex i);
  double computeMean(FloatTensorVector& vector);
};

} /* namespace combigrid */
} /* namespace sgpp */
