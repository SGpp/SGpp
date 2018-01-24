// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>

#include <sgpp/combigrid/algebraic/FirstMomentNormStrategy.hpp>
#include <sgpp/combigrid/algebraic/VarianceNormStrategy.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class PolynomialStochasticCollocation : public CombigridSurrogateModel {
 public:
  PolynomialStochasticCollocation(sgpp::combigrid::CombigridSurrogateModelConfiguration& config);
  virtual ~PolynomialStochasticCollocation();

  double mean() override;
  double variance() override;

  void getComponentSobolIndices(sgpp::base::DataVector& componentSsobolIndices,
                                bool normalized = true) override;
  void getTotalSobolIndices(sgpp::base::DataVector& totalSobolIndices,
                            bool normalized = true) override;

  void updateOperation(
      std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation) override;
  void updateOperation(
      std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridOperation) override;
  void updateOperation(
      std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridOperation) override;
  void updateConfig(sgpp::combigrid::CombigridSurrogateModelConfiguration config) override;

 private:
  void initializeTensorOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::shared_ptr<LevelManager> levelManager);

  void initializeBounds();
  void initializeWeightFunctions();
  void initializeNormStrategies();

  bool updateStatus();
  double computeMean();
  double computeVariance();

  // global polynomial basis
  std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> legendreBasis;
  // orthogonal basis for pdf values
  sgpp::combigrid::WeightFunctionsCollection weightFunctions;

  size_t numGridPoints;
  sgpp::combigrid::FloatTensorVector expansionCoefficients;

  // mean and variance storage
  bool computedMeanFlag;
  double ev;
  std::unique_ptr<sgpp::combigrid::FirstMomentNormStrategy> firstMomentNormstrategy;

  bool computedVarianceFlag;
  double var;
  std::unique_ptr<sgpp::combigrid::VarianceNormStrategy> varianceNormStrategy;
};

} /* namespace combigrid */
} /* namespace sgpp */
