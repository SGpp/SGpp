// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/algebraic/NormStrategy.hpp>
#include <sgpp/combigrid/grid/distribution/LejaPointDistribution.hpp>
#include <sgpp/combigrid/grid/distribution/UniformBoundaryPointDistribution.hpp>
#include <sgpp/combigrid/grid/growth/LinearGrowthStrategy.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/hierarchy/NonNestedPointHierarchy.hpp>
#include <sgpp/combigrid/grid/ordering/IdentityPointOrdering.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/LevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridSummationStrategy.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <memory>
#include <vector>

#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluationStrategy.hpp>

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
 * For method documentation, please refer to CombigridMultiOperation.
 */
class CombigridOperation {
  // unique_ptr would be possible, but gives SWIG errors
  std::shared_ptr<CombigridOperationImpl> impl;

 public:
  CombigridOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, MultiFunction func, bool exploitNesting = true,
      FullGridSummationStrategyType summationStrategyType = FullGridSummationStrategyType::LINEAR,
      std::shared_ptr<NormStrategy<FloatScalarVector>> normStrategy = nullptr);

  CombigridOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage,
      FullGridSummationStrategyType summationStrategyType = FullGridSummationStrategyType::LINEAR,
      std::shared_ptr<NormStrategy<FloatScalarVector>> normStrategy = nullptr);

  CombigridOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, GridFunction gridFunc, bool exploitNesting = true,
      FullGridSummationStrategyType summationStrategyType = FullGridSummationStrategyType::LINEAR,
      std::shared_ptr<NormStrategy<FloatScalarVector>> normStrategy = nullptr);

  void setParameters(base::DataVector const &param = base::DataVector(0));  // clears automatically

  double getResult();

  double evaluate(size_t q, base::DataVector const &param = base::DataVector(0));
  std::shared_ptr<LevelManager> getLevelManager();
  void setLevelManager(std::shared_ptr<LevelManager> levelManager);

  /**
   * @return the storage containing the computed function values at evaluation points.
   */
  std::shared_ptr<AbstractCombigridStorage> getStorage();

  /**
   *  @param storage the storage containing the coefficients precalculated by some other operation
   */
  void setStorage(std::shared_ptr<AbstractCombigridStorage> storage);

  /**
   * @return the point hierarchies containing the grid points in each direction
   */
  std::vector<std::shared_ptr<AbstractPointHierarchy>> getPointHierarchies();

  /**
   * @return the evaluator prototypes in each direction
   */
  std::shared_ptr<AbstractFullGridEvaluator<FloatScalarVector>> getFullGridEval();

  /**
   * @return the number of function values that have been computed via this CombigridOperation
   * during its lifetime. For a nested grid, this number matches numGridPoints() if only one
   * computation is performed, i.e. no previous data has been cleared via evaluate() or
   * setParameters(). Its computation is not optimized, but currently faster than
   * numGridPoints().
   */
  size_t numStoredFunctionValues();

  /**
   * @return the total number of different (multi-dimensional) grid points that have been used for
   * the current evaluation. This number is reset when clear() is called on the CombigridEvaluator
   * via evaluate() or setParameters().
   * This method is currently not optimized and can be slow!
   */
  size_t numGridPoints();

  /**
   * @return the number of dimensions
   */
  size_t numDims();

  /**
   * @return An upper bound for the number of points (function evaluations) used for the current
   * computation. This bound is exact if nesting is used or if otherwise each grid point only occurs
   * in exactly one level.
   */
  size_t getUpperPointBound() const;

  static std::shared_ptr<CombigridOperation> createExpClenshawCurtisPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createExpChebyshevPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createExpLejaPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createExpL2LejaPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createExpUniformPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createExpUniformBoundaryPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createExpUniformLinearInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createExpUniformBoundaryLinearInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createLinearClenshawCurtisPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createLinearLejaPolynomialInterpolation(
      size_t numDimensions, MultiFunction func, size_t growthFactor = 2);
  static std::shared_ptr<CombigridOperation> createLinearL2LejaPolynomialInterpolation(
      size_t numDimensions, MultiFunction func, size_t growthFactor = 2);
  static std::shared_ptr<CombigridOperation> createLinearUniformPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createLinearUniformBoundaryPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);
  static std::shared_ptr<CombigridOperation> createLinearLejaQuadrature(size_t numDimensions,
                                                                        MultiFunction func,
                                                                        size_t growthFactor = 2);
  static std::shared_ptr<CombigridOperation> createLinearL2LejaQuadrature(size_t numDimensions,
                                                                          MultiFunction func,
                                                                          size_t growthFactor = 2);
  static std::shared_ptr<CombigridOperation> createExpClenshawCurtisQuadrature(size_t numDimensions,
                                                                               MultiFunction func);

  static std::shared_ptr<CombigridOperation> auxiliaryBsplineFunction(
      size_t numDimensions, MultiFunction func, sgpp::combigrid::CombiHierarchies::Collection grids,
      sgpp::combigrid::CombiEvaluators::Collection evaluators, size_t degree);

  static std::shared_ptr<CombigridOperation> createExpUniformBoundaryBsplineInterpolation(
      size_t numDimensions, MultiFunction func, size_t degree);
  static std::shared_ptr<CombigridOperation> createExpClenshawCurtisBsplineInterpolation(
      size_t numDimensions, MultiFunction func, size_t degree);
  static std::shared_ptr<CombigridOperation> createExpChebyshevBsplineInterpolation(
      size_t numDimensions, MultiFunction func, size_t degree);
  static std::shared_ptr<CombigridOperation> createLinearLejaBsplineInterpolation(
      size_t numDimensions, MultiFunction func, size_t degree, size_t growthFactor);
  static std::shared_ptr<CombigridOperation> createLinearL2LejaBsplineInterpolation(
      size_t numDimensions, MultiFunction func, size_t degree, size_t growthFactor);
  static std::shared_ptr<CombigridOperation> createExpUniformBoundaryBsplineQuadrature(
      size_t numDimensions, MultiFunction func, size_t degree);
};
} /* namespace combigrid */
} /* namespace sgpp*/
