// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/FloatArrayVector.hpp>
#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>

#include <memory>

namespace sgpp {
namespace combigrid {

class CombiHierarchies {
 public:
  static std::shared_ptr<AbstractPointHierarchy> linearLeja(size_t growthFactor = 2);
  static std::shared_ptr<AbstractPointHierarchy> expLeja();
  static std::shared_ptr<AbstractPointHierarchy> expUniform();
  static std::shared_ptr<AbstractPointHierarchy> expClenshawCurtis();
};

class CombiEvaluators {
 public:
  static std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> polynomialInterpolation();
  static std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> linearInterpolation();
  static std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> quadrature();
  static std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> multiPolynomialInterpolation();
  static std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> multiLinearInterpolation();
  static std::shared_ptr<AbstractLinearEvaluator<FloatArrayVector>> multiQuadrature();
};

} /* namespace combigrid */
} /* namespace sgpp */
