// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>

namespace sgpp {
namespace base {

/**
 * Abstract class representing a system of linear equations.
 * All row and column indices are zero based.
 */
class SLE {
 public:
  /**
   * Constructor.
   */
  SLE() {}

  /**
   * Destructor.
   */
  virtual ~SLE() {}

  /**
   * Pure virtual method for checking if a matrix entry vanishes or not.
   *
   * @param i     row index
   * @param j     column index
   * @return      whether the (i,j)-th entry of the matrix is non-zero
   */
  virtual bool isMatrixEntryNonZero(size_t i, size_t j) = 0;

  /**
   * Pure virtual method for retrieving a matrix entry.
   *
   * @param i     row index
   * @param j     column index
   * @return      (i,j)-th entry of the matrix
   */
  virtual double getMatrixEntry(size_t i, size_t j) = 0;

  /**
   * Multiply the matrix with a vector.
   * Standard implementation with \f$\mathcal{O}(n^2)\f$ scalar
   * multiplications.
   *
   * @param       x   vector to be multiplied
   * @param[out]  y   \f$y = Ax\f$
   */
  virtual void matrixVectorMultiplication(const base::DataVector& x, base::DataVector& y) {
    const size_t n = getDimension();
    y.resize(n);
    y.setAll(0.0);

    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        y[i] += getMatrixEntry(i, j) * x[j];
      }
    }
  }

  /**
   * Count all non-zero entries.
   * Standard implementation with \f$\mathcal{O}(n^2)\f$ checks.
   *
   * @return number of non-zero entries
   */
  virtual size_t countNNZ() {
    const size_t n = getDimension();
    size_t nnz = 0;

    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        if (isMatrixEntryNonZero(i, j)) {
          nnz++;
        }
      }
    }

    return nnz;
  }

  /**
   * Pure virtual method returning the dimension (number of rows/columns)
   * of the system.
   *
   * @return  system dimension
   */
  virtual size_t getDimension() const = 0;

  /**
   * @return whether this system derives from CloneableSLE or not
   *         (standard: false)
   */
  virtual bool isCloneable() const { return false; }
};
}  // namespace base
}  // namespace sgpp
