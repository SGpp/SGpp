// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONQUADRATURE_HPP
#define OPERATIONQUADRATURE_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * This class provides the quadrature of a sparse grid function
 */
class OperationQuadrature {
 public:
  /**
   * Constructor
   */
  OperationQuadrature() {}

  /**
   * Destructor
   */
  virtual ~OperationQuadrature() {}

  /**
   * Integrate the sparse grid function
   *
   * @param alpha the function's values in the nodal basis
   */
  virtual double doQuadrature(DataVector& alpha) = 0;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONQUADRATURE_HPP */
