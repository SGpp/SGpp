// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMATRIX_HPP
#define OPERATIONMATRIX_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * Abstract definition of a matrix operator interface.
 * Every time you need to apply a matrix to the ansatzfunction's
 * coefficients derive a class from OperationMatrix
 */
class OperationMatrix {
 public:
  /**
   * Constructor
   */
  OperationMatrix() {}

  /**
   * Destructor
   */
  virtual ~OperationMatrix() {}

  /**
   * starts the Multiplication with the matrix
   *
   * @param alpha DataVector that contains the ansatzfunctions' coefficients
   * @param result DataVector into which the result of the Laplace operation is stored
   */
  virtual void mult(DataVector& alpha, DataVector& result) = 0;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONMATRIX_HPP */
