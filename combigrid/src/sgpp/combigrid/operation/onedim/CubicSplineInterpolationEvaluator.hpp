// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class CubicSplineInterpolationEvaluator : public AbstractLinearEvaluator<FloatScalarVector> {
 public:
  CubicSplineInterpolationEvaluator();
  ~CubicSplineInterpolationEvaluator() override;
  CubicSplineInterpolationEvaluator(CubicSplineInterpolationEvaluator const &other);

  std::vector<FloatScalarVector> getBasisValues() override { return basisValues; }
  std::vector<double> getBasisCoefficients() override { return basisCoefficients; }

  void setGridPoints(std::vector<double> const &newXValues) override;
  void setBasisCoefficientsAtGridPoints(std::vector<double> &functionValues) override;
  std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> cloneLinear() override;
  bool needsOrderedPoints() override;
  bool needsParameter() override;
  void setParameter(FloatScalarVector const &param) override;

 private:
  void computeBasisValues();

  double evaluationPoint;
  std::vector<FloatScalarVector> basisValues;
  std::vector<double> basisCoefficients;
  std::vector<double> xValues;

  struct SplineSet {
    double a;
    double b;
    double c;
    double d;
  };

  std::vector<std::vector<SplineSet>> gridCoefficients;
};

} /* namespace combigrid */
} /* namespace sgpp */
