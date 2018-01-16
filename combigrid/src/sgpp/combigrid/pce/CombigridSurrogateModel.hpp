// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>

namespace sgpp {
namespace combigrid {

enum class CombigridSurrogateModelsType {
  POLYNOMIAL_CHAOS_EXPANSION,
  POLYNOMIAL_STOCHASTIC_COLLOCATION,
  BSPLINE_STOCHASTIC_COLLOCATION
};

class CombigridSurrogateModelConfiguration {
 public:
  // type
  CombigridSurrogateModelsType type;

  // operation
  std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation;
  std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation;
  std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation;

  // structure
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies;
  std::shared_ptr<AbstractCombigridStorage> storage;
  std::shared_ptr<LevelManager> levelManager;
  //  std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> levelStructure;

  // basis function for tensor operation
  std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction;
  sgpp::combigrid::OrthogonalBasisFunctionsCollection basisFunctions;

  // bounds for stochastic collocation
  sgpp::base::DataVector bounds;

  // Bspline degree
  size_t degree;

  //  Bspline coefficients
  std::shared_ptr<AbstractCombigridStorage> coefficientStorage;

  size_t numDims;
};

// --------------------------------------------------------------------------

class CombigridSurrogateModel {
 public:
  CombigridSurrogateModel(sgpp::combigrid::CombigridSurrogateModelConfiguration& config);
  virtual ~CombigridSurrogateModel();

  virtual double mean() = 0;
  virtual double variance() = 0;
  virtual void getComponentSobolIndices(sgpp::base::DataVector& componentSsobolIndices,
                                        bool normalized = true) = 0;
  virtual void getTotalSobolIndices(sgpp::base::DataVector& totalSobolIndices,
                                    bool normalized = true) = 0;

  virtual void updateOperation(
      std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation) = 0;
  virtual void updateOperation(
      std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridOperation) = 0;
  virtual void updateOperation(
      std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridOperation) = 0;

  virtual void updateConfig(sgpp::combigrid::CombigridSurrogateModelConfiguration config) = 0;

  sgpp::combigrid::CombigridSurrogateModelConfiguration& getConfig();

 protected:
  sgpp::combigrid::CombigridSurrogateModelConfiguration config;

  size_t numDims;
};

} /* namespace combigrid */
} /* namespace sgpp */
