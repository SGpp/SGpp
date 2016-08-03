// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/MultiFunction.hpp>
#include <sgpp/combigrid/algebraic/FloatArrayVector.hpp>
#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/multidim/LevelManager.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

class CombigridMultiOperationImpl;
// we use pimpl for not having to include all the template stuff
// in the header

/**
 * Interface class for simple usage of the combigrid module. Via a CombigridMultiOperation, the
 * evaluation (at multiple interpolation points at once) of the computation pipeline can be easily
 * managed.
 * There are two main ways to create this class:
 * - The point hierarchies and evaluators etc. are created by the user and passed to the
 * constructor
 * - One of the static methods is used. They provide some sensible isotropic configurations.
 *
 * Via the LevelManager, which can be get and set, one can control which adaptivity criterion might
 * be used. For easy evaluation, there is an evaluate()-method, which does all the work at once and
 * generates a regular level structure. To get more control over the level structure, one may
 * proceed as follows:
 * - Set the parameters for interpolation via setParameters()
 * - add combigrid levels via getLevelManager()->some_add_levels_function()
 * - fetch the result via getResult().
 */
class CombigridMultiOperation {
  std::shared_ptr<CombigridMultiOperationImpl> impl;  // unique_ptr causes SWIG errors

 public:
  CombigridMultiOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, MultiFunction func);

  CombigridMultiOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager,
      std::shared_ptr<AbstractCombigridStorage> storage);

  // TODO(holzmudd): add extra functions, for example for configuring the storage

  /**
   * Sets the parameters for upcoming computations and clears the data structures (removes old
   * computed data)
   */
  void setParameters(std::vector<base::DataVector> const &params);

  /**
   * Sets the parameters for upcoming computations and clears the data structures (removes old
   * computed data). The evaluation points are assumed to be the columns of the matrix.
   */
  void setParameters(base::DataMatrix const &params);

  base::DataVector getResult();

  base::DataVector evaluate(size_t q, std::vector<base::DataVector> const &params);

  std::shared_ptr<LevelManager> getLevelManager();
  void setLevelManager(std::shared_ptr<LevelManager> levelManager);

  std::shared_ptr<AbstractMultiStorage<FloatArrayVector>> getDifferences();

  // TODO(holzmudd): add static constructor functions
  static std::shared_ptr<CombigridMultiOperation> createExpClenshawCurtisPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridMultiOperation> createLinearClenshawCurtisPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridMultiOperation> createExpLejaPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridMultiOperation> createLinearLejaPolynomialInterpolation(
      size_t numDimensions, MultiFunction func, size_t growthFactor = 2);
  static std::shared_ptr<CombigridMultiOperation> createLinearLejaQuadrature(
      size_t numDimensions, MultiFunction func, size_t growthFactor = 2);
  static std::shared_ptr<CombigridMultiOperation> createExpUniformLinearInterpolation(
      size_t numDimensions, MultiFunction func);
};
} /* namespace combigrid */
} /* namespace sgpp*/
