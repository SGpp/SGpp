// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVALGRADIENT_HPP
#define OPERATIONNAIVEEVALGRADIENT_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <vector>

namespace sgpp {
namespace base {

/**
 * Abstract operation for evaluating a linear combination of basis functions and its gradient.
 * The "naive" is indicating that classes implementing this operation should use a "naive"
 * approach, e.g. by evaluating all basis functions by brute force.
 */
class OperationNaiveEvalGradient {
 public:
  /**
   * Constructor.
   */
  OperationNaiveEvalGradient() {
  }

  /**
   * Destructor.
   */
  virtual ~OperationNaiveEvalGradient() {
  }

  /**
   * @param       alpha     coefficient vector
   * @param       point     evaluation point
   * @param[out]  gradient  gradient vector of the linear combination
   * @return                value of the linear combination
   */
  virtual double evalGradient(const DataVector& alpha,
                               const DataVector& point,
                               DataVector& gradient) = 0;

  /**
   * @param       alpha     coefficient matrix (each column is a coefficient vector)
   * @param       point     evaluation point
   * @param[out]  value     values of the linear combination
   * @param[out]  gradient  Jacobian of the linear combination (each row is a gradient vector)
   */
  virtual void evalGradient(const DataMatrix& alpha,
                            const DataVector& point,
                            DataVector& value,
                            DataMatrix& gradient) {
    const size_t d = point.getSize();
    const size_t m = alpha.getNcols();
    DataVector curAlpha(alpha.getNrows());
    DataVector curGradient(d);

    value.resize(m);
    gradient.resize(m, d);

    for (size_t j = 0; j < m; j++) {
      alpha.getColumn(j, curAlpha);
      value[j] = evalGradient(curAlpha, point, curGradient);
      gradient.setRow(j, curGradient);
    }
  }
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONNAIVEEVALGRADIENT_HPP */
