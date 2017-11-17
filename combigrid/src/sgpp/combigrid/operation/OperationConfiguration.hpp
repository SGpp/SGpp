// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>

namespace sgpp {
namespace combigrid {

enum class CombiEvaluatorTypes {
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

  // interpolation
  Tensor_PolynomialInterpolation
};

struct OperationConfiguration {
  CombiEvaluatorTypes type;
  size_t degree = 3;
  std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis;

  OperationConfiguration(CombiEvaluatorTypes type) : type(type){};
  OperationConfiguration(CombiEvaluatorTypes type, size_t degree) : type(type), degree(degree){};
  OperationConfiguration(CombiEvaluatorTypes type, size_t degree,
                         std::shared_ptr<AbstractInfiniteFunctionBasis1D> functionBasis)
      : type(type), degree(degree), functionBasis(functionBasis){};

  ~OperationConfiguration(){};
};

} /* namespace combigrid */
} /* namespace sgpp*/
