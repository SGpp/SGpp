// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONREGULARIZATIONDIAGONAL_HPP
#define OPERATIONREGULARIZATIONDIAGONAL_HPP

#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

/**
 * Implementation of the application of a diagonal matrix to a
 * DataVector for regularization.
 * This class implements several scaling possibilities.
 */
class OperationRegularizationDiagonal : public base::OperationMatrix {
 protected:
  // to remember mode
  int mode;
  // to remember parameter k
  float_t k;
  // to remember state of grid in terms of number of grid points
  size_t size;
  // ro remember grid's storage
  base::GridStorage* storage;
  // to store diagonal
  base::DataVector diagonal;

  // Diagonal entries have been provided
  // static const int DIAGMATRIX = 0;

  /**
   * Check for mode and run corresponding init
   */
  virtual void init();

  /**
   * Initialize Hkmix.
   * Diagonal entries are
   * @f$\prod_{k=1}^d \langle \phi_k(x_k),\phi_k(x_k) \rangle_{H^k}@f$.
   * @param k Parameter k
   */
  virtual void initHkmix(float_t k) = 0;

  /**
   * Initialize H0HkLaplace.
   * Diagonal entries are
   * @f$\sum_{k=1}^d \langle \phi_k(x_k),\phi_k(x_k) \rangle_{H^k}
   * \prod_{l\neq k} \langle \phi_k(x_k),\phi_k(x_k) \rangle_{H^0}@f$.
   * @param k Parameter k
   */
  virtual void initH0HkLaplace(float_t k) = 0;

  /**
   * Initialize ISOTROPIC_PENALTY, ignores constructor parameter k.
   * Diagonal entries are
   * @f$\frac{1}{\max\{l_1,\dots,l_d\}-\min\{l_1,\dots,l_d\}+1}d@f$.
   */
  virtual void initIsotropicPenalty();

  /**
   * Initialize ANISOTROPIC_PENALTY, ignores constructor parameter k.
   * Diagonal entries are
   * @f$\frac{1}{2}\log(1+(\frac{\max\{l_1,\dots,l_d\}}{\max\{\min\{l_1,\dots,l_d\},1\}})d)@f$.
   */
  virtual void initAnisotropicPenalty();

 public:
  /// Diagonal scaling for @f$H^k_\textbf{mix}@f$-norm (product of
  /// @f$H^k@f$ 1D norms)
  static const int HKMIX = 1;
  /// Diagonal scaling for pseudo-laplace, replacing @f$L_2@f$ with
  /// @f$H^0@f$ in 1d
  static const int H0HKLAPLACE = 2;
  // Penalizes strongly isotropic subspaces more
  static const int ISOTROPIC_PENALTY = 3;
  // Penalizes strongly anisotropic subspaces more
  static const int ANISOTROPIC_PENALTY = 4;

  /**
   * Constructor of OperationRegularizationDiagonal.
   * Constructor should most likely call init() in subclasses.
   * @param storage Pointer to grid's storage object
   * @param mode Mode, specifying which regularization to use. Example:
   * OperationRegularizationDiagonal::HKMIX.
   * @param k Parameter for @f$H^k@f$
   */
  OperationRegularizationDiagonal(base::GridStorage* storage, int mode, float_t k);

  /**
   * Destructor
   */
  virtual ~OperationRegularizationDiagonal();

  /**
   * Multiplication with diagonal matrix.
   * @param alpha Source vector
   * @param result Result vector (entries are overwritten)
   */
  virtual void mult(base::DataVector& alpha, base::DataVector& result);
};
}  // namespace datadriven
}  // namespace SGPP
#endif /* OPERATIONREGULARIZATIONDIAGONAL_HPP */
