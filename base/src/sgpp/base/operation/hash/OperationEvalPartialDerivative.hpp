// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVALPARTIALDERIVATIVE_HPP
#define OPERATIONEVALPARTIALDERIVATIVE_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <vector>

namespace sgpp {
namespace base {

/**
 * Abstract operation for evaluating a partial derivative of a linear combination of
 * basis functions.
 * The "naive" is indicating that classes implementing this operation should use a "naive"
 * approach, e.g. by evaluating all basis functions by brute force.
 */
class OperationEvalPartialDerivative {
 public:
  /**
   * Constructor.
   */
  OperationEvalPartialDerivative() {
  }

  /**
   * Destructor.
   */
  virtual ~OperationEvalPartialDerivative() {
  }

  /**
   * @param       alpha               coefficient vector
   * @param       point               evaluation point
   * @param       derivDim            dimension in which the partial derivative should be taken
   *                                  (0, ..., d-1)
   * @param[out]  partialDerivative   value of the partial derivative of the linear combination
   * @return                          value of the linear combination
   */
  virtual double evalPartialDerivative(const DataVector& alpha,
                                        const DataVector& point,
                                        size_t derivDim,
                                        double& partialDerivative) = 0;

  /**
   * @param       alpha               coefficient matrix (each column is a coefficient vector)
   * @param       point               evaluation point
   * @param       derivDim            dimension in which the partial derivative should be taken
   *                                  (0, ..., d-1)
   * @param[out]  value               values of the linear combination
   * @param[out]  partialDerivative   values of the partial derivatives of the linear combination
   *                                  (the j-th entry corresponds to the j-th column of alpha)
   */
  virtual void evalPartialDerivative(const DataMatrix& alpha,
                                     const DataVector& point,
                                     size_t derivDim,
                                     DataVector& value,
                                     DataVector& partialDerivative) {
    const size_t m = alpha.getNcols();
    DataVector curAlpha(alpha.getNrows());

    value.resize(m);
    partialDerivative.resize(m);

    for (size_t j = 0; j < m; j++) {
      alpha.getColumn(j, curAlpha);
      value[j] = evalPartialDerivative(curAlpha, point, derivDim, partialDerivative[j]);
    }
  }
  /// untransformed evaluation point (temporary vector)
  DataVector pointInUnitCube;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONEVALPARTIALDERIVATIVE_HPP */
