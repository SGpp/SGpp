// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef QUADRATUREEVALUATOR_HPP_
#define QUADRATUREEVALUATOR_HPP_

#include <sgpp/combigrid/operation/onedim/AbstractLinearEvaluator.hpp>
#include <sgpp/combigrid/algebraic/ScalarVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/SingleFunction.hpp>

#include <vector>
#include <functional>

namespace sgpp {
namespace combigrid {

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
