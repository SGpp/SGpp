// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONDOTPRODUCTLINEAR_HPP
#define OPERATIONDOTPRODUCTLINEAR_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

class OperationDotProductLinear {
 public:
  /**
   * Constructor of OperationTestLinear
   *
   * @param storage Pointer to the grid's gridstorage obejct
   */
  explicit OperationDotProductLinear(sgpp::base::GridStorage* storage) : storage(storage) {}

  /**
   * Destructor
   */
  virtual ~OperationDotProductLinear() {}

  virtual double eval(sgpp::base::DataVector& x1, sgpp::base::DataVector& x2);

 protected:
  /// Pointer to the grid's GridStorage object
  sgpp::base::GridStorage* storage;
};
}  // namespace datadriven
}  // namespace sgpp

#endif /* OPERATIONDOTPRODUCTLINEAR_HPP */
