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

class FirstMomentNormStrategy : public NormStrategy<FloatTensorVector> {
 public:
  FirstMomentNormStrategy(
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction, size_t numDims,
      sgpp::combigrid::SingleFunction weightFunction, bool isOrthogonal,
      sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  FirstMomentNormStrategy(
      std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction,
      sgpp::combigrid::WeightFunctionsCollection& weightFunctions, bool isOrthogonal,
      sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  FirstMomentNormStrategy(sgpp::combigrid::OrthogonalBasisFunctionsCollection& basisFunctions,
                          sgpp::combigrid::WeightFunctionsCollection& weightFunctions,
                          bool isOrthogonal,
                          sgpp::base::DataVector const& bounds = sgpp::base::DataVector(0));

  ~FirstMomentNormStrategy() override;

  /**
   * computes the first moment of the given tensor vector
   *
   * @param vector tensor
   * @return first moment
   */
  double norm(FloatTensorVector& vector) override;

 private:
  bool isOrthogonal;
  sgpp::base::DataVector bounds;

  sgpp::combigrid::OrthogonalBasisFunctionsCollection basisFunctions;
  sgpp::combigrid::WeightFunctionsCollection weightFunctions;

  std::map<MultiIndex, double> lookupTable;

  double quad(MultiIndex i, GaussLegendreQuadrature& quadRule);
  double computeMean(FloatTensorVector& vector);

  void initializeBounds();
};

} /* namespace combigrid */
} /* namespace sgpp */
