// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVALHESSIAN_HPP
#define OPERATIONNAIVEEVALHESSIAN_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <vector>

namespace SGPP {
namespace base {

/**
 * Abstract operation for evaluating a linear combination of basis functions, its gradient
 * and its Hessian.
 * The "naive" is indicating that classes implementing this operation should use a "naive"
 * approach, e.g. by evaluating all basis functions by brute force.
 */
class OperationNaiveEvalHessian {
 public:
  /**
   * Constructor.
   */
  OperationNaiveEvalHessian() {
  }

  /**
   * Destructor.
   */
  virtual ~OperationNaiveEvalHessian() {
  }

  /**
   * Pure virtual method for evaluating a linear combination of basis functions, its gradient
   * and its Hessian.
   *
   * @param       alpha       coefficient vector
   * @param       point       evaluation point
   * @param[out]  gradient    gradient vector of the linear combination
   * @param[out]  hessian     Hessian matrix of the linear combination
   * @return                  value of the linear combination
   */
  float_t evalHessian(const DataVector& alpha,
                      const std::vector<float_t>& point,
                      DataVector& gradient,
                      DataMatrix& hessian) {
    DataVector p(point);
    return evalHessian(alpha, p, gradient, hessian);
  }

  /**
   * Convenience function for using DataVector as points.
   *
   * @param       alpha       coefficient vector
   * @param       point       evaluation point
   * @param[out]  gradient    gradient vector of the linear combination
   * @param[out]  hessian     Hessian matrix of the linear combination
   * @return                  value of the linear combination
   */
  virtual float_t evalHessian(const DataVector& alpha,
                              const DataVector& point,
                              DataVector& gradient,
                              DataMatrix& hessian) = 0;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONNAIVEEVALHESSIAN_HPP */
