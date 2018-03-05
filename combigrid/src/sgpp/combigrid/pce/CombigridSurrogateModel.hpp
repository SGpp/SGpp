// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/functions/ProbabilityDensityFunction1D.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

enum class CombigridSurrogateModelsType {
  POLYNOMIAL_CHAOS_EXPANSION,
  POLYNOMIAL_STOCHASTIC_COLLOCATION,
  BSPLINE_STOCHASTIC_COLLOCATION
};

class CombigridSurrogateModelConfiguration {
 public:
  CombigridSurrogateModelConfiguration();
  virtual ~CombigridSurrogateModelConfiguration();

  // type
  CombigridSurrogateModelsType type;

  // structure
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies;
  std::shared_ptr<AbstractCombigridStorage> storage;
  std::shared_ptr<sgpp::combigrid::LevelManager> levelManager;
  std::shared_ptr<sgpp::combigrid::TreeStorage<uint8_t>> levelStructure;

  // basis function for tensor operation
  std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction;
  sgpp::combigrid::OrthogonalBasisFunctionsCollection basisFunctions;
  std::shared_ptr<sgpp::combigrid::FloatTensorVector> expansionCoefficients;

  // extract knowledge directly from operation
  std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> tensorOperation;

  // bounds for stochastic collocation
  sgpp::base::DataVector bounds;

  // Bspline degree
  size_t degree;
  //  Bspline coefficients
  std::shared_ptr<AbstractCombigridStorage> coefficientStorage;

  // ToDo (rehmemk) initialize with a single weight function and build a homogene
  // WeightFunctionsCollection from that. Just like basisFunction above
  // weight functions
  sgpp::combigrid::WeightFunctionsCollection weightFunctions;
  std::shared_ptr<sgpp::combigrid::SingleFunction> weightFunction;

  bool enableLevelManagerStatsCollection;
  size_t numDimensions;

  void loadFromCombigridOperation(std::shared_ptr<CombigridOperation> op,
                                  bool loadLevelStructure = true);
  void loadFromCombigridOperation(std::shared_ptr<CombigridMultiOperation> op,
                                  bool loadLevelStructure = true);
  void loadFromCombigridOperation(std::shared_ptr<CombigridTensorOperation> op,
                                  bool loadLevelStructure = true);
};

// --------------------------------------------------------------------------

class CombigridSurrogateModel {
 public:
  explicit CombigridSurrogateModel(sgpp::combigrid::CombigridSurrogateModelConfiguration& config);
  virtual ~CombigridSurrogateModel();

  virtual double eval(sgpp::base::DataVector& x) = 0;
  virtual void eval(sgpp::base::DataMatrix& xs, sgpp::base::DataVector& res) = 0;

  virtual double mean() = 0;
  virtual double variance() = 0;
  virtual void getComponentSobolIndices(sgpp::base::DataVector& componentSsobolIndices,
                                        bool normalized = true) = 0;
  virtual void getTotalSobolIndices(sgpp::base::DataVector& totalSobolIndices,
                                    bool normalized = true) = 0;

  virtual void updateConfig(sgpp::combigrid::CombigridSurrogateModelConfiguration config) = 0;

  virtual size_t numGridPoints() = 0;
  virtual std::shared_ptr<LevelInfos> getInfoOnAddedLevels() = 0;

  sgpp::combigrid::CombigridSurrogateModelConfiguration& getConfig();

 protected:
  sgpp::combigrid::CombigridSurrogateModelConfiguration config;

  size_t numDims;
};

} /* namespace combigrid */
} /* namespace sgpp */
