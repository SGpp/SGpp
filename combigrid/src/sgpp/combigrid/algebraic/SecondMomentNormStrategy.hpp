// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/NormStrategy.hpp>
#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>

#include <vector>
#include <map>

namespace sgpp {
namespace combigrid {

class SecondMomentNormStrategy : public NormStrategy<FloatTensorVector> {
 public:
  SecondMomentNormStrategy(
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction, size_t numDims,
      sgpp::combigrid::SingleFunction weightFunction, bool isOrthogonal,
      sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  SecondMomentNormStrategy(
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction,
      sgpp::combigrid::WeightFunctionsCollection& weightFunctions, bool isOrthogonal,
      sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  SecondMomentNormStrategy(sgpp::combigrid::OrthogonalBasisFunctionsCollection& basisFunctions,
                           sgpp::combigrid::WeightFunctionsCollection& weightFunctions,
                           bool isOrthogonal,
                           sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  ~SecondMomentNormStrategy() override;

  /**
   * computes the second moment of the given tensor vector
   *
   * @param vector tensor
   * @return second moment
   */
  double norm(FloatTensorVector& vector) override;

 private:
  bool isOrthogonal;
  sgpp::base::DataVector bounds;
  std::map<MultiIndex, double> innerProducts;

  sgpp::combigrid::OrthogonalBasisFunctionsCollection basisFunctions;
  sgpp::combigrid::WeightFunctionsCollection weightFunctions;

  double quad(sgpp::combigrid::MultiIndex i, sgpp::combigrid::MultiIndex j,
              GaussLegendreQuadrature& quadRule);
  double computeSecondMoment(sgpp::combigrid::FloatTensorVector& vector);

  void initializeBounds();

  void joinMultiIndices(MultiIndex& ix, MultiIndex& jx, MultiIndex& kx);
};

} /* namespace combigrid */
} /* namespace sgpp */
