// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/MultiFunction.hpp>
#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/UniformPointDistribution.hpp>
#include <sgpp/combigrid/grid/growth/LinearGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NonNestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/IdentityPointOrdering.hpp>
#include <sgpp/combigrid/operation/multidim/LevelManager.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

class CombigridOperationImpl;
// we use pimpl for not having to include all the template stuff in
// the header

/**
 * Interface class for simple usage of the combigrid module. Via a CombigridOperation, the
 * evaluation (interpolation at a single point or quadrature) of the computation pipeline can be
 * easily managed.
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
 *
 * For method documentation, refer to CombigridMultiOperation.
 */
class CombigridOperation {
  std::shared_ptr<CombigridOperationImpl>
      impl;  // unique_ptr would be possible, but gives SWIG errors

 public:
  CombigridOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, MultiFunction func);

  CombigridOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager,
      std::shared_ptr<AbstractCombigridStorage> storage);

  void setParameters(base::DataVector const &param = base::DataVector(0));  // clears automatically

  double getResult();

  double evaluate(size_t q, base::DataVector const &param = base::DataVector(0));
  std::shared_ptr<LevelManager> getLevelManager();
  void setLevelManager(std::shared_ptr<LevelManager> levelManager);

  std::shared_ptr<AbstractCombigridStorage> getStorage();

  /**
   * @return the number of function values that have been computed via this CombigridOperation
   * during its lifetime. For a nested grid, this number matches numGridPoints() if only one
   * computation is performed, i.e. no previous data has been cleared via evaluate() or
   * setParameters(). Its computation is not optimized, but currently faster than numGridPoints().
   */
  size_t numStoredFunctionValues();

  /**
   * @return the total number of different (multi-dimensional) grid points that have been used for
   * the current evaluation. This number is reset when clear() is called on the CombigridEvaluator
   * via evaluate() or setParameters().
   * This method is currently not optimized and can be slow!
   */
  size_t numGridPoints();

  static std::shared_ptr<CombigridOperation> createExpClenshawCurtisPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createExpChebyshevPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createExpLejaPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createExpUniformPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createLinearClenshawCurtisPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createLinearLejaPolynomialInterpolation(
      size_t numDimensions, MultiFunction func, size_t growthFactor = 2);
  static std::shared_ptr<CombigridOperation> createLinearUniformPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createLinearLejaQuadrature(size_t numDimensions,
                                                                        MultiFunction func,
                                                                        size_t growthFactor = 2);
  static std::shared_ptr<CombigridOperation> createExpClenshawCurtisQuadrature(size_t numDimensions,
                                                                               MultiFunction func);
  static std::shared_ptr<CombigridOperation> createExpUniformLinearInterpolation(
      size_t numDimensions, MultiFunction func);
};
} /* namespace combigrid */
} /* namespace sgpp*/
