// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/tools/sle/system/CloneableSLE.hpp>

#include <cstddef>

namespace sgpp {
namespace base {

/**
 * Full linear system, essentially a wrapper around base::DataMatrix.
 */
class FullSLE : public CloneableSLE {
 public:
  /**
   * Constructor.
   * Do not destruct the matrix A before this object!
   *
   * @param A     coefficient matrix
   */
  explicit FullSLE(base::DataMatrix& A) : CloneableSLE(), A(A) {}

  /**
   * Destructor.
   */
  ~FullSLE() override {}

  /**
   * @param i     row index
   * @param j     column index
   * @return      whether the (i,j)-th entry of the matrix is non-zero
   */
  inline bool isMatrixEntryNonZero(size_t i, size_t j) override { return (A(i, j) != 0.0); }

  /**
   * @param i     row index
   * @param j     column index
   * @return      (i,j)-th entry of the matrix
   */
  inline double getMatrixEntry(size_t i, size_t j) override { return A(i, j); }

  /**
   * @return  coefficient matrix
   */
  base::DataMatrix& getA() { return A; }

  size_t getDimension() const override { return A.getNrows(); }

  /**
   * Clones the linear system.
   * Because A is stored as a reference, A is not copied (only b).
   *
   * @param[out] clone pointer to cloned object
   */
  void clone(std::unique_ptr<CloneableSLE>& clone) const override {
    clone = std::unique_ptr<CloneableSLE>(new FullSLE(A));
  }

 protected:
  /// coefficient matrix
  base::DataMatrix& A;
};
}  // namespace base
}  // namespace sgpp
