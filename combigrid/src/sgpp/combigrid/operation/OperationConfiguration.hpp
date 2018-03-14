// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/GeneralFunction.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>

#include <iostream>
#include <memory>

namespace sgpp {
namespace combigrid {

enum class CombiEvaluatorTypes {
  // dummy for construction of evaluators without configuration.
  NO_TYPE,

  // interpolation
  Scalar_PolynomialInterpolation,
  Scalar_LinearInterpolation,
  Scalar_CubicSplineInterpolation,
  Scalar_BSplineInterpolation,

  Multi_PolynomialInterpolation,
  Multi_LinearInterpolation,
  Multi_CubicSplineInterpolation,
  Multi_BSplineInterpolation,

  // quadrature
  Scalar_PolynomialQuadrature,
  Scalar_BSplineQuadrature,

  Multi_PolynomialQuadrature,
  Multi_BSplineQuadrature,

  // scalar product
  Multi_BSplineScalarProduct,
  Multi_PolynomialScalarProduct,

  // interpolation
  Tensor_PolynomialInterpolation,
  Tensor_BSplineInterpolation
};

struct EvaluatorConfiguration {
  CombiEvaluatorTypes type;
  size_t degree = 3;
  std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis;

  EvaluatorConfiguration() : EvaluatorConfiguration(CombiEvaluatorTypes::NO_TYPE, 0, nullptr) {}
  EvaluatorConfiguration(CombiEvaluatorTypes type, size_t degree,
                         std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis)
      : type(type), degree(degree), functionBasis(functionBasis) {}

  explicit EvaluatorConfiguration(CombiEvaluatorTypes type)
      : EvaluatorConfiguration(type, 0, nullptr) {}
  EvaluatorConfiguration(CombiEvaluatorTypes type, size_t degree)
      : EvaluatorConfiguration(type, degree, nullptr) {}
  EvaluatorConfiguration(CombiEvaluatorTypes type,
                         std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis)
      : EvaluatorConfiguration(type, 0, functionBasis) {}

  ~EvaluatorConfiguration() {}
};

} /* namespace combigrid */
} /* namespace sgpp*/
