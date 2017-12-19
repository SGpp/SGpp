// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/NormStrategy.hpp>
#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>

namespace sgpp {
namespace combigrid {

class SecondMomentNormStrategy : public NormStrategy<FloatTensorVector> {
 public:
  SecondMomentNormStrategy(
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction, size_t numDims,
      bool isOrthogonal);
  SecondMomentNormStrategy(
      std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& basisFunctions,
      bool isOrthogonal);
  virtual ~SecondMomentNormStrategy();

  double norm(FloatTensorVector& vector) {
    if (isOrthogonal) {
      double sum = 0.0;
      for (auto it = vector.getValues()->getStoredDataIterator(); it->isValid(); it->moveToNext()) {
        double coeff = it->value().value();
        sum += coeff * coeff;
      }
      return sum;
    } else {
      return computeSecondMoment(vector);
    }
  }

 private:
  bool isOrthogonal;
  std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>> basisFunctions;

  double quad(sgpp::combigrid::MultiIndex i, sgpp::combigrid::MultiIndex j);
  double computeSecondMoment(sgpp::combigrid::FloatTensorVector& vector);
};

} /* namespace combigrid */
} /* namespace sgpp */
