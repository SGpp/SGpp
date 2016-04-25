// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONCONVERT_HPP
#define OPERATIONCONVERT_HPP

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>


namespace sgpp {
namespace base {

/**
 * Operation that cpnverts a given basis into the normal, linear hat basis and vice versa
 *
 */
class OperationConvert {
 public:
  /**
   * Constructor
   */
  OperationConvert() {}

  /**
   * Destructor
   */
  virtual ~OperationConvert() {}

  /**
   * Convert given basis into linear hat basis.
   */
  virtual void doConvertToLinear(DataVector& alpha) = 0;


  /**
   * Convert from a linear coefficient vector into given basis.
   */
  virtual void doConvertFromLinear(DataVector& alpha) = 0;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONCONVERT_HPP */
