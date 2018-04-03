// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/algebraic/FloatArrayVector.hpp>
#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/algebraic/FloatTensorVector.hpp>
#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/OperationConfiguration.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This class provides standard configurations of point hierarchies. The methods names allude to the
 * used growth strategy and point distribution. All hierarchies provided here are nested.
 */
class CombiHierarchies {
 public:
  static std::shared_ptr<AbstractPointHierarchy> linearLeja(size_t growthFactor = 2);
  static std::shared_ptr<AbstractPointHierarchy> linearLeja(SingleFunction weightFunction,
                                                            size_t growthFactor = 2);
  static std::shared_ptr<AbstractPointHierarchy> linearL2Leja(size_t growthFactor = 2);
  static std::shared_ptr<AbstractPointHierarchy> linearL2Leja(SingleFunction weightFunction,
                                                              size_t growthFactor = 2,
                                                              size_t numAdditionalPoints = 10);
  static std::shared_ptr<AbstractPointHierarchy> expLeja();
  static std::shared_ptr<AbstractPointHierarchy> expLeja(SingleFunction weightFunction);
  static std::shared_ptr<AbstractPointHierarchy> expL2Leja();
  static std::shared_ptr<AbstractPointHierarchy> expL2Leja(SingleFunction weightFunction,
                                                           size_t numAdditionalPoints = 10);
  static std::shared_ptr<AbstractPointHierarchy> expUniform();
  static std::shared_ptr<AbstractPointHierarchy> expClenshawCurtis();
  static std::shared_ptr<AbstractPointHierarchy> expChebyshev();
  static std::shared_ptr<AbstractPointHierarchy> expUniformBoundary();

  /**
   * Not efficient because it is not nested.
   */
  static std::shared_ptr<AbstractPointHierarchy> linearUniform(size_t growthFactor = 2);

  /**
   * Not efficient because it is not nested.
   */
  static std::shared_ptr<AbstractPointHierarchy> linearClenshawCurtis(size_t growthFactor = 2);
  /**
   * Not efficient because it is not nested.
   */
  static std::shared_ptr<AbstractPointHierarchy> linearChebyshev(size_t growthFactor = 2);

  /**
   * Not efficient because it is not nested.
   */
  static std::shared_ptr<AbstractPointHierarchy> linearUniformBoundary(size_t growthFactor = 2);

  typedef std::vector<std::shared_ptr<AbstractPointHierarchy>> Collection;
};

/**
 * This class provides standard configurations of 1D-evaluators (single- and multi-evaluation).
 */
class CombiEvaluators {
 public:
  // scalar evaluators
  static std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> createCombiScalarEvaluator(
      EvaluatorConfiguration operationConfig);

  static std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> polynomialInterpolation();
  static std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> linearInterpolation();
  static std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> cubicSplineInterpolation();
  static std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> BSplineInterpolation(
      size_t degree);
  static std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> quadrature();

  static std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> BSplineQuadrature(
      size_t degree);
  static std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> BSplineQuadrature(
      size_t degree, sgpp::combigrid::SingleFunction weight_function, size_t numAdditionalPoints,
      bool normalizeWeights);
  static std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> BSplineQuadrature(
      size_t degree, sgpp::combigrid::SingleFunction weight_function, size_t numAdditionalPoints,
      double a, double b, bool normalizeWeights);

  // array evaluators
  static std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> createCombiMultiEvaluator(
      EvaluatorConfiguration operationConfig);

  static std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> BSplineScalarProduct(
      size_t degree);
  static std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> polynomialScalarProduct();

  static std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> multiPolynomialInterpolation();
  static std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> multiLinearInterpolation();
  static std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> multiCubicSplineInterpolation();
  static std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> multiBSplineInterpolation(
      size_t degree);
  static std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> multiBSplineQuadrature(
      size_t degree);
  static std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> multiQuadrature();
  static std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> multiQuadrature(
      SingleFunction func, bool normalizeWeights);

  // tensor evaluators
  static std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>> createCombiTensorEvaluator(
      EvaluatorConfiguration operationConfig);

  static std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>> tensorInterpolation(
      std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis);

  static std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>> tensorBSplineInterpolation(
      size_t degree = 3);

  typedef std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>>> Collection;
  typedef std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>>> MultiCollection;
  typedef std::vector<std::shared_ptr<AbstractLinearEvaluator<FloatTensorVector>>> TensorCollection;
};

} /* namespace combigrid */
} /* namespace sgpp */
