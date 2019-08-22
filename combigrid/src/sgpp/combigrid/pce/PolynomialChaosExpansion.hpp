// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#ifdef USE_DAKOTA
#include <OrthogPolyApproximation.hpp>
#endif

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class PolynomialChaosExpansion : public CombigridSurrogateModel {
 public:
  explicit PolynomialChaosExpansion(sgpp::combigrid::CombigridSurrogateModelConfiguration& config);
  ~PolynomialChaosExpansion() override;

  double eval(sgpp::base::DataVector& x) override;
  void eval(sgpp::base::DataMatrix& xs, sgpp::base::DataVector& res) override;

  double mean() override;
  double variance() override;

  void getComponentSobolIndices(sgpp::base::DataVector& componentSsobolIndices,
                                bool normalized = true) override;
  void getTotalSobolIndices(sgpp::base::DataVector& totalSobolIndices,
                            bool normalized = true) override;

  void updateConfig(sgpp::combigrid::CombigridSurrogateModelConfiguration config) override;

  size_t numGridPoints() override;
  std::shared_ptr<LevelInfos> getInfoOnAddedLevels() override;

 private:
  bool updateStatus();
  void computeComponentSobolIndices();

#ifdef USE_DAKOTA
  std::shared_ptr<Pecos::OrthogPolyApproximation> orthogPoly;
#endif

  // basis functions
  sgpp::combigrid::OrthogonalBasisFunctionsCollection basisFunctions;

  // tensor operation
  std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation;
  sgpp::combigrid::FloatTensorVector expansionCoefficients;

  size_t currentNumGridPoints;
  bool computedSobolIndicesFlag;
  sgpp::base::DataVector sobolIndices;
};

} /* namespace combigrid */
} /* namespace sgpp */
