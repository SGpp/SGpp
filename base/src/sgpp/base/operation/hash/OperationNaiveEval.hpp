// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONNAIVEEVAL_HPP
#define OPERATIONNAIVEEVAL_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>


namespace sgpp {
namespace base {

/**
 * Operation that evaluates the current sparse grid function defined
 * by the coefficient vector @em alpha at a given point.
 * In contrast to OperationEval, implementations of OperationNaiveEval should
 * use a "naive" method for evaluating sparse grid functions, e.g. evaluate
 * all basis functions by brute force.
 */
class OperationNaiveEval {
 public:
  /**
   * Default constructor.
   */
  OperationNaiveEval() {}

  /**
   * Destructor
   */
  virtual ~OperationNaiveEval() {}

  /**
   * @param alpha     coefficient vector
   * @param point     evaluation point
   * @return          value of the linear combination
   */
  virtual double eval(const DataVector& alpha,
                      const DataVector& point) = 0;

  /**
   * @param      alpha  coefficient matrix (each column is a coefficient vector)
   * @param      point  evaluation point
   * @param[out] value  values of the linear combination
   */
  virtual void eval(const DataMatrix& alpha,
                    const DataVector& point,
                    DataVector& value) {
    const size_t m = alpha.getNcols();
    DataVector curAlpha(alpha.getNrows());

    value.resize(m);

    for (size_t j = 0; j < m; j++) {
      alpha.getColumn(j, curAlpha);
      value[j] = eval(curAlpha, point);
    }
  }
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONNAIVEEVAL_HPP */
