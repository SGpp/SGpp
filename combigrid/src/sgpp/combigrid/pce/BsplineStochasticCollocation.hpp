// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/operation/multidim/sparsegrid/LTwoScalarProductHashMapNakBsplineBoundaryCombigrid.hpp>
#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class BsplineStochasticCollocation : public CombigridSurrogateModel {
 public:
  explicit BsplineStochasticCollocation(
      sgpp::combigrid::CombigridSurrogateModelConfiguration& config);
  ~BsplineStochasticCollocation() override;

  void eval(sgpp::base::DataMatrix& xs, sgpp::base::DataVector& res) override;
  double eval(sgpp::base::DataVector& x) override;

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
  void initializeOperations(std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
                            std::shared_ptr<AbstractCombigridStorage> storage,
                            std::shared_ptr<LevelManager> levelManager);

  bool updateStatus();
  double computeMean();
  double computeVariance();

  void countPolynomialTerms();
  size_t additionalQuadraturePoints(OrthogonalPolynomialBasisType polyType);

  double quad(sgpp::combigrid::MultiIndex i);
  double quad(sgpp::combigrid::MultiIndex i, sgpp::combigrid::MultiIndex j);

  // pdf values
  sgpp::combigrid::WeightFunctionsCollection weightFunctions;
  bool customWeightFunction;

  size_t currentNumGridPoints;

  // mean and variance storage
  bool computedMeanFlag;
  double ev;
  bool computedVarianceFlag;
  double var;
  // basis coefficients for Bspline interpolation
  std::shared_ptr<sgpp::combigrid::AbstractCombigridStorage> coefficientStorage;
  //  ToDo(rehmemk)
  size_t numthreads = 4;
  LTwoScalarProductHashMapNakBsplineBoundaryCombigrid scalarProducts;

  // operation
  std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation;
  std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation;
};

} /* namespace combigrid */
} /* namespace sgpp */
