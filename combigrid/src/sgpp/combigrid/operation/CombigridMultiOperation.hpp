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
  /**
   * Constructs a CombigridMultiOperation with the given hierarchies, evaluators, level manager and
   * function to evaluate.
   */
  CombigridMultiOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager, MultiFunction func);

  /**
   * Constructs a CombigridMultiOperation with the given hierarchies, evaluators, level manager and
   * CombigridStorage that contains the function to evaluate.
   */
  CombigridMultiOperation(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies,
      std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> evaluatorPrototypes,
      std::shared_ptr<LevelManager> levelManager,
      std::shared_ptr<AbstractCombigridStorage> storage);

  /**
   * Sets the parameters for upcoming computations and clears the data structures (removes old
   * computed data). This is only relevant for methods with parameters (e.g. interpolation, but not
   * quadrature).
   * @param params The parameters at which the function should be evaluated.
   */
  void setParameters(std::vector<base::DataVector> const &params);

  /**
   * Sets the parameters for upcoming computations and clears the data structures (removes old
   * computed data). This is only relevant for methods with parameters (e.g. interpolation, but not
   * quadrature). The evaluation points are assumed to be the columns of the matrix.
   */
  void setParameters(base::DataMatrix const &params);

  /**
   * Returns the current computed value.
   */
  base::DataVector getResult();

  /**
   * Computes the result with regular levels up to 1-norm q (levels start from zero) and with
   * parameter params. This is a convenience function. If you need other functionality, use
   * getLevelManager() and operate directly on the LevelManager.
   */
  base::DataVector evaluate(size_t q, std::vector<base::DataVector> const &params);

  /**
   * See the other version of evaluate() and the setParameters() overloads.
   */
  base::DataVector evaluate(size_t q, base::DataMatrix const &params);

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
  std::shared_ptr<AbstractMultiStorage<FloatArrayVector>> getDifferences();

  /**
   * @return the number of function values that have been computed via this CombigridMultiOperation
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
   * Returns a CombigridMultiOperation doing polynomial interpolation on a Clenshaw-Curtis grid with
   * an exponential growth (nested points).
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be interpolated.
   */
  static std::shared_ptr<CombigridMultiOperation> createExpClenshawCurtisPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);

  /**
   * Returns a CombigridMultiOperation doing polynomial interpolation on a Chebyshev grid with
   * an exponential growth (nested points).
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be interpolated.
   */
  static std::shared_ptr<CombigridMultiOperation> createExpChebyshevPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);

  /**
   * Returns a CombigridMultiOperation doing polynomial interpolation on a Clenshaw-Curtis grid with
   * a linear growth (not nested points).
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be interpolated.
   */
  static std::shared_ptr<CombigridMultiOperation> createLinearClenshawCurtisPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);

  /**
   * Returns a CombigridMultiOperation doing polynomial interpolation on a Leja grid with
   * an exponential growth (nested points).
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be interpolated.
   */
  static std::shared_ptr<CombigridMultiOperation> createExpLejaPolynomialInterpolation(
      size_t numDimensions, MultiFunction func);

  /**
   * Returns a CombigridMultiOperation doing polynomial interpolation on a Leja grid with
   * linear growth (nested points).
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be interpolated.
   * @param growthFactor Parameter for the linear growth strategy. For level l, 1 + growthFactor * l
   * points are used.
   */
  static std::shared_ptr<CombigridMultiOperation> createLinearLejaPolynomialInterpolation(
      size_t numDimensions, MultiFunction func, size_t growthFactor = 2);

  /**
   * Returns a CombigridMultiOperation doing quadrature (based on integrals of Lagrange polynomials)
   * on a Leja grid with linear growth (nested points).
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be integrated.
   * @param growthFactor Parameter for the linear growth strategy. For level l, 1 + growthFactor * l
   * points are used.
   */
  static std::shared_ptr<CombigridMultiOperation> createLinearLejaQuadrature(
      size_t numDimensions, MultiFunction func, size_t growthFactor = 2);

  /**
   * Returns a CombigridMultiOperation doing quadrature (based on integrals of Lagrange polynomials)
   * on a Clenshaw-Curtis with exponential growth (nested points).
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be integrated.
   */
  static std::shared_ptr<CombigridMultiOperation> createExpClenshawCurtisQuadrature(
      size_t numDimensions, MultiFunction func);

  /**
   * Returns a CombigridMultiOperation doing (multi-)linear interpolation
   * on a uniform grid with boundary using an exponential growth strategy (nested points).
   * @param numDimensions Dimensionality of the problem.
   * @param func Function to be interpolated.
   */
  static std::shared_ptr<CombigridMultiOperation> createExpUniformLinearInterpolation(
      size_t numDimensions, MultiFunction func);
};
} /* namespace combigrid */
} /* namespace sgpp*/
