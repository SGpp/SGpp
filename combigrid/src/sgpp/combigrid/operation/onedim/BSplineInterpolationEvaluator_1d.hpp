// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>


#include <vector>

namespace sgpp {
namespace combigrid {

class BSplineInterpolationEvaluator_1d : public AbstractLinearEvaluator<FloatScalarVector> {
 public:
  BSplineInterpolationEvaluator_1d();
  BSplineInterpolationEvaluator_1d(size_t degree);
  virtual ~BSplineInterpolationEvaluator_1d();
  BSplineInterpolationEvaluator_1d(BSplineInterpolationEvaluator_1d const &other);

  std::vector<FloatScalarVector> getBasisValues() override { return basisValues; }
  std::vector<double> getBasisCoefficients() override { return basisCoefficients; }

  void setDegree(size_t const &deg);
  void setGridPoints(std::vector<double> const &newXValues) override;
  void setFunctionValuesAtGridPoints(std::vector<double> &functionValues) override;
  std::shared_ptr<AbstractLinearEvaluator<FloatScalarVector>> cloneLinear() override;
  bool needsOrderedPoints() override;
  bool needsParameter() override;
  void setParameter(FloatScalarVector const &param) override;

 private:
  void computeBasisValues();
  //void computeCoefficients();

  double evaluationPoint;
  std::vector<FloatScalarVector> basisValues;
  std::vector<double> basisCoefficients;
  std::vector<double> xValues;

  size_t degree;
};

} /* namespace combigrid */
} /* namespace sgpp */
