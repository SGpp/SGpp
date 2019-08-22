// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_POLYNOMIALINTERPOLATIONEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_POLYNOMIALINTERPOLATIONEVALUATOR_HPP_

#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This evaluator does polynomial interpolation (using the barycentric approach) on the given grid
 * points.
 */
class PolynomialInterpolationEvaluator : public AbstractLinearEvaluator<FloatScalarVector> {
  double evaluationPoint;
  std::vector<FloatScalarVector> basisValues;
  std::vector<double> basisCoefficients;

  std::vector<double> wValues;
  std::vector<double> xValues;

  void computeBasisValues();

 public:
  PolynomialInterpolationEvaluator();
  ~PolynomialInterpolationEvaluator() override;
  PolynomialInterpolationEvaluator(PolynomialInterpolationEvaluator const &other);

  std::vector<FloatScalarVector> getBasisValues() override { return basisValues; }
  std::vector<double> getBasisCoefficients() override { return basisCoefficients; }

  void setGridPoints(std::vector<double> const &newXValues) override;
  void setBasisCoefficientsAtGridPoints(std::vector<double> &functionValues) override;
  std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> cloneLinear() override;
  bool needsOrderedPoints() override;
  bool needsParameter() override;
  void setParameter(FloatScalarVector const &param) override;

  // TODO(holzmudd): eval could be optimized...
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_ONEDIM_POLYNOMIALINTERPOLATIONEVALUATOR_HPP_ */
