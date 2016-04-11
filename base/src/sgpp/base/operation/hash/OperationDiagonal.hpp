// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONDIAGONAL_HPP
#define OPERATIONDIAGONAL_HPP

#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Implementation of identity Operation for all kinds of grids
 */
class OperationDiagonal : public OperationMatrix {
 private:
  DataVector diag;

 public:
  /**
   * Constructor of OperationDiagonal
   */
  explicit OperationDiagonal(DataVector v) : diag(v) {}

  /**
   * Destructor
   */
  ~OperationDiagonal() override {}

  void mult(DataVector& alpha, DataVector& result) override {
    result = DataVector(alpha);  // componentwise_mult isn't a pure func.
    result.componentwise_mult(diag);
  }
};

}  // namespace base
}  // namespace sgpp
#endif /* OPERATIONDIAGONAL_HPP */
