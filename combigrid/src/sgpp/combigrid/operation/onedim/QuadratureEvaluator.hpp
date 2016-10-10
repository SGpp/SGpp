// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef QUADRATUREEVALUATOR_HPP_
#define QUADRATUREEVALUATOR_HPP_

#include <sgpp/combigrid/SingleFunction.hpp>
#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>

#include <functional>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This evaluator does quadrature based on the given grid points. The quadrature weights are
 * obtained by (numerically) integrating the Lagrange polynomials on the given grid points.
 * In the constructor, a weight function may be passed whose values at the grid points are
 * multiplied with the given function values.
 */
class QuadratureEvaluator : public AbstractLinearEvaluator<FloatScalarVector> {
  std::vector<double> xValues;
  std::vector<FloatScalarVector> weights;
  sgpp::combigrid::SingleFunction weight_function;
  bool normalizeWeights;

 public:
  QuadratureEvaluator(
      sgpp::combigrid::SingleFunction weight_function =
          sgpp::combigrid::SingleFunction(constantFunction<double>(static_cast<double>(1.0))),
      bool normalizeWeights = false);
  QuadratureEvaluator(QuadratureEvaluator const &other);
  virtual ~QuadratureEvaluator();

  virtual std::vector<FloatScalarVector> getBasisCoefficients() { return weights; }

  virtual void setGridPoints(std::vector<double> const &newXValues);

  virtual bool needsOrderedPoints();
  virtual bool needsParameter();
  virtual void setParameter(FloatScalarVector const &param);

  virtual std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> cloneLinear();
};

#endif /* QUADRATUREEVALUATOR_HPP_ */

} /* namespace combigrid */
} /* namespace sgpp*/
