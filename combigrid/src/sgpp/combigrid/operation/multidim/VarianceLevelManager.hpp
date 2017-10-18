// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/operation/multidim/LevelManager.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluator.hpp>
#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class VarianceLevelManager : public LevelManager {
 public:
  explicit VarianceLevelManager(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::vector<OrthogonalPolynomialBasisType> basisFunctionsType);
  virtual ~VarianceLevelManager();

  virtual std::shared_ptr<LevelManager> clone();

 protected:
  double computePriority(MultiIndex const &level) override;

 private:
  /**
   * Point hierarchies that provide the grid points in each dimension.
   */
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies;

  /**
   * Function values storage.
   */
  std::shared_ptr<AbstractCombigridStorage> storage;

  std::shared_ptr<AbstractFullGridEvaluator<FloatTensorVector>> fullGridEval;
};

} /* namespace combigrid */
} /* namespace sgpp */
