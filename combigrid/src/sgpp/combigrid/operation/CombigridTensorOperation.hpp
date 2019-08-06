// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/algebraic/NormStrategy.hpp>
#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/multidim/LevelManager.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>

#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <memory>
#include <vector>

#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractFullGridEvaluationStrategy.hpp>

namespace sgpp {
namespace combigrid {

class CombigridTensorOperationImpl;

class CombigridTensorOperation {
  std::shared_ptr<CombigridTensorOperationImpl> impl;  // unique_ptr causes SWIG errors

 public:
  /**
   * Constructs a CombigridTensorOperation with the given hierarchies, evaluators, level manager
   * and
   * function to evaluate.
   */
  CombigridTensorOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, MultiFunction func, bool exploitNesting = true,
      FullGridSummationStrategyType summationStrategyType = FullGridSummationStrategyType::LINEAR,
      std::shared_ptr<NormStrategy<FloatTensorVector>> normStrategy = nullptr);

  /**
   * Constructs a CombigridTensorOperation with the given hierarchies, evaluators, level manager and
   * CombigridStorage that contains the function to evaluate.
   */
  CombigridTensorOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, std::shared_ptr<AbstractCombigridStorage> storage,
      FullGridSummationStrategyType summationStrategyType = FullGridSummationStrategyType::LINEAR,
      std::shared_ptr<NormStrategy<FloatTensorVector>> normStrategy = nullptr);

  /**
   * Constructs a CombigridTensorOperation with the given hierarchies, evaluators, level manager
   * and grid function together with the information whether the same point on different levels
   * should be able to have different function values.
   */
  CombigridTensorOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, GridFunction gridFunc, bool exploitNesting,
      FullGridSummationStrategyType summationStrategyType = FullGridSummationStrategyType::LINEAR,
      std::shared_ptr<NormStrategy<FloatTensorVector>> normStrategy = nullptr);

  /**
   * Sets the parameters for upcoming computations and clears the data structures (removes old
   * computed data). This is only relevant for methods with parameters.
   * @param params The parameters at which the function should be evaluated.
   */
  void setParameters(std::vector<FloatTensorVector> const &params);

  /**
   * Returns the current computed value.
   */
  FloatTensorVector getResult();

  /**
   * Computes the result with regular levels up to 1-norm q (levels start from zero) and with
   * parameter params. This is a convenience function. If you need other functionality, use
   * getLevelManager() and operate directly on the LevelManager.
   */
  FloatTensorVector evaluate(
      size_t q, std::vector<FloatTensorVector> const &params = std::vector<FloatTensorVector>());

  /**
   * @return the storage containing the computed function values at evaluation points.
   */
  std::shared_ptr<AbstractCombigridStorage> getStorage();

  /**
   * Via the LevelManager, more options are available than are provided directly by this class.
   */
  std::shared_ptr<LevelManager> getLevelManager();

  /**
   * Can be used to set the level manager, e.g. if one of the static constructor functions has been
   * used.
   */
  void setLevelManager(std::shared_ptr<LevelManager> levelManager);

  /**
   * Returns a storage of the differences (Deltas) that have been computed by the
   * CombigridEvaluator.
   */
  std::shared_ptr<AbstractMultiStorage<FloatTensorVector>> getDifferences();

  /**
   * @return the number of function values that have been computed via this CombigridTensorOperation
   * during its lifetime. For a nested grid, this number matches numGridPoints() if only one
   * computation is  performed, i.e. no previous data has been cleared via evaluate() or
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

  /**
   * @return the point hierarchies containing the grid points in each direction
   */
  std::vector<std::shared_ptr<AbstractPointHierarchy>> getPointHierarchies();

  /**
   * Returns a CombigridTensorOperation doing polynomial interpolation on a Clenshaw-Curtis grid
   * with
   * an exponential growth (nested points).
   * @param functionBasis global basis to which the result should be transformed
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be interpolated.
   */
  static std::shared_ptr<CombigridTensorOperation> createExpClenshawCurtisPolynomialInterpolation(
      std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
      MultiFunction func);

  /**
   * Returns a CombigridTensorOperation doing polynomial interpolation on a Chebyshev grid with
   * an exponential growth (nested points).
   * @param functionBasis global basis to which the result should be transformed
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be interpolated.
   */
  static std::shared_ptr<CombigridTensorOperation> createExpChebyshevPolynomialInterpolation(
      std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
      MultiFunction func);

  /**
   * Returns a CombigridTensorOperation doing polynomial interpolation on a Clenshaw-Curtis grid
   * with a linear growth (not nested points).
   * @param functionBasis global basis to which the result should be transformed
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be interpolated.
   */
  static std::shared_ptr<CombigridTensorOperation>
  createLinearClenshawCurtisPolynomialInterpolation(
      std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
      MultiFunction func);

  /**
   * Returns a CombigridTensorOperation doing polynomial interpolation on a Leja grid with
   * an exponential growth (nested points).
   * @param functionBasis global basis to which the result should be transformed
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be interpolated.
   */
  static std::shared_ptr<CombigridTensorOperation> createExpLejaPolynomialInterpolation(
      std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
      MultiFunction func);

  /**
   * Returns a CombigridTensorOperation doing polynomial interpolation on a L2Leja grid with
   * an exponential growth (nested points).
   * @param functionBasis global basis to which the result should be transformed
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be interpolated.
   */
  static std::shared_ptr<CombigridTensorOperation> createExpL2LejaPolynomialInterpolation(
      std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
      MultiFunction func);

  /**
   * Returns a CombigridTensorOperation doing polynomial interpolation on a Leja grid with
   * linear growth (nested points).
   * @param functionBasis global basis to which the result should be transformed
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be interpolated.
   * @param growthFactor Parameter for the linear growth strategy. For level l, 1 + growthFactor * l
   * points are used.
   */
  static std::shared_ptr<CombigridTensorOperation> createLinearLejaPolynomialInterpolation(
      std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
      MultiFunction func, size_t growthFactor = 2);

  /**
   * Returns a CombigridTensorOperation doing polynomial interpolation on a L2Leja grid with
   * linear growth (nested points).
   * @param functionBasis global basis to which the result should be transformed
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be interpolated.
   * @param growthFactor Parameter for the linear growth strategy. For level l, 1 + growthFactor * l
   * points are used.
   */
  static std::shared_ptr<CombigridTensorOperation> createLinearL2LejaPolynomialInterpolation(
      std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis, size_t numDimensions,
      MultiFunction func, size_t growthFactor = 2);

  /**
   * Transforms the basic structures of an arbitrary operation to a tensor operation
   *
   * @param pointHierarchies univariate grids
   * @param storage function value storage
   * @param functionBasis global basis function to which the result should be transformed
   * @param summationStrategyType strategy to gather the results of the univariate evaluators on
   * @return tensor operation with the same grid as given by the parameters
   */
  static std::shared_ptr<CombigridTensorOperation> createOperationTensorPolynomialInterpolation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::shared_ptr<AbstractCombigridStorage> storage,
      std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis,
      FullGridSummationStrategyType summationStrategyType = FullGridSummationStrategyType::LINEAR);

  /**
   * Transforms the basic structures of an arbitrary operation to a tensor operation
   *
   * @param pointHierarchies univariate grids
   * @param storage function value storage
   * @param functionBases vector of global basis functions to which the result should be transformed
   * @param summationStrategyType strategy to gather the results of the univariate evaluators on
   * @return tensor operation with the same grid as given by the parameters
   */
  static std::shared_ptr<CombigridTensorOperation> createOperationTensorPolynomialInterpolation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::shared_ptr<AbstractCombigridStorage> storage,
      OrthogonalBasisFunctionsCollection &functionBases,
      FullGridSummationStrategyType summationStrategyType = FullGridSummationStrategyType::LINEAR);

  /**
   * Transforms the basic structures of an arbitrary operation to a tensor operation
   *
   * @param pointHierarchies univariate grids
   * @param storage function value storage
   * @param levelManager provides level structures that are copied to the new tensor operation
   * @param summationStrategyType strategy to gather the results of the univariate evaluators on
   * @return tensor operation with the same grid as given by the parameters
   */
  static std::shared_ptr<CombigridTensorOperation> createOperationTensorBSplineInterpolation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::shared_ptr<AbstractCombigridStorage> storage, std::shared_ptr<LevelManager> levelManager,
      FullGridSummationStrategyType summationStrategyType = FullGridSummationStrategyType::LINEAR);

  static std::shared_ptr<CombigridTensorOperation> createExpUniformBoundaryBSplineInterpolation(
      size_t numDimensions, MultiFunction func, size_t degree);
};

} /* namespace combigrid */
} /* namespace sgpp */
