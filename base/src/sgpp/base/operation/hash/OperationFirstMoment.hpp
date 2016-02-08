// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONFIRSTMOMENT_HPP
#define OPERATIONFIRSTMOMENT_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * This class provides the first moment of a sparse grid function
 */
class OperationFirstMoment {
 public:
  /**
   * Constructor
   */
  OperationFirstMoment() {}

  /**
   * Destructor
   */
  virtual ~OperationFirstMoment() {}

  /**
   * Integrate the sparse grid function
   *
   * @param alpha the function's values in the nodal basis
   */
  virtual float_t doQuadrature(const DataVector& alpha) = 0;

};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONFIRSTMOMENT_HPP */
