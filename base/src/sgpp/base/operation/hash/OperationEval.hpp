// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONEVAL_HPP
#define OPERATIONEVAL_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>


namespace sgpp {
namespace base {

/**
 * Operation that evaluates the current sparse grid function defined
 * by the coefficient vector @em alpha at a given point.
 */
class OperationEval {
 public:
  /**
   * Default constructor.
   */
  OperationEval() {}

  /**
   * Destructor
   */
  virtual ~OperationEval() {}

  /**
   * Evaluates the sparse grid function at a given point.
   *
   * @param alpha The coefficients of the sparse grid's basis functions
   * @param point The coordinates of the evaluation point
   */
  double eval(const DataVector& alpha,
               const std::vector<double>& point) {
    DataVector p(point);
    return eval(alpha, p);
  }

  /**
   * Evaluates the sparse grid function at a given point.
   *
   * @param alpha The coefficients of the sparse grid's basis functions
   * @param point The coordinates of the evaluation point
   */
  virtual double eval(const DataVector& alpha,
                       const DataVector& point) = 0;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONEVAL_HPP */
