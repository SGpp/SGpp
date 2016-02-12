// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONGAMMALOGLINEARBOUNDARY_HPP
#define OPERATIONGAMMALOGLINEARBOUNDARY_HPP

#include <sgpp/pde/algorithm/UpDownTwoOpDims.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace finance {

/**
 * Implements the Gamma-Operation, that is needed
 * the solve the multidimensional Black Scholes
 * (in the non-transformed case, this corresponds to matrix E in Master's thesis)
 * equation
 *
 */
class OperationGammaLogLinearBoundary : public SGPP::pde::UpDownTwoOpDims {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's SGPP::base::GridStorage object
   * @param coef vector that contains the constant coefficients of this operation
   */
  OperationGammaLogLinearBoundary(SGPP::base::GridStorage* storage, SGPP::base::DataMatrix& coef);

  /**
   * Destructor
   */
  virtual ~OperationGammaLogLinearBoundary();

 protected:
  /**
   * Up-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
   * Applies the up-part of the one-dimensional mass matrix in one dimension.
   * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i < l_j} \alpha_j \phi_j(x) dx.\f]
   *
   * @param dim dimension in which to apply the up-part
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  virtual void up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

  /**
   * Down-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
   * Applies the down-part of the one-dimensional mass matrix in one dimension.
   * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i\geq l_j} \alpha_j \phi_j(x) dx.\f]
   *
   * @param dim dimension in which to apply the down-part
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  virtual void down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim);

  /**
   * down-Gradient step in dimension <i>dim</i> applies the x phi dphi operation
   * in one dimension
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that down-Gradient is applied
   */
  virtual void downOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                            size_t dim);

  /**
   * up-Gradient step in dimension <i>dim</i> applies the x phi dphi operation
   * in one dimension
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that up-Gradient is applied
   */
  virtual void upOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                          size_t dim);

  /**
   * down-Gradient step in dimension <i>dim</i> applies the x dphi phi operation
   * in one dimension
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that down-Gradient is applied
   */
  virtual void downOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                            size_t dim);

  /**
   * up-Gradient step in dimension <i>dim</i> applies the x dphi phi operation
   * in one dimension
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that up-Gradient is applied
   */
  virtual void upOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                          size_t dim);

  /**
   * down-Gradient multiplied with a squared x step in dimension <i>dim</i> applies the x^2 dphi
   * dphi operation
   * in one dimension
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that down-Gradient is applied
   */
  void downOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                               size_t dim);

  /**
   * up-Gradient multiplied with a squared x step in dimension <i>dim</i> applies the x^2 dphi dphi
   * operation
   * in one dimension
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that up-Gradient is applied
   */
  void upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                             size_t dim);
};
}  // namespace finance
}  // namespace SGPP

#endif /* OPERATIONGAMMALOGLINEARBOUNDARY_HPP */
