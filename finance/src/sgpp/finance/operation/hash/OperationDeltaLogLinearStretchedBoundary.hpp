// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONDELTALOGLINEARSTRETCHEDBOUNDARY_HPP
#define OPERATIONDELTALOGLINEARSTRETCHEDBOUNDARY_HPP

#include <sgpp/pde/algorithm/UpDownOneOpDim.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace finance {

/**
 * Implements the Delta-Operation, that is needed
 * the solve the multidimensional Black Scholes
 * (in the non-transformed case, this corresponds to matrix F in Master's thesis)
 * equation
 *
 */
class OperationDeltaLogLinearStretchedBoundary : public sgpp::pde::UpDownOneOpDim {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's sgpp::base::GridStorage object
   * @param coef reference to a sgpp::base::DataVector object that contains the bilinear form's
   * constant coefficients
   */
  OperationDeltaLogLinearStretchedBoundary(sgpp::base::GridStorage* storage,
                                           sgpp::base::DataVector& coef);

  /**
   * Destructor
   */
  virtual ~OperationDeltaLogLinearStretchedBoundary();

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
  virtual void up(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  /**
   * Down-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
   * Applies the down-part of the one-dimensional mass matrix in one dimension.
   * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i\geq l_j} \alpha_j \phi_j(x) dx.\f]
   *
   * @param dim dimension in which to apply the down-part
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  virtual void down(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  /**
   * down-Gradient step in dimension <i>dim</i> applies the x dphi phi operation
   * in one dimension
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that down-Gradient is applied
   */
  virtual void downOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);

  /**
   * up-Gradient step in dimension <i>dim</i> applies the x dphi phi operation
   * in one dimension
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that up-Gradient is applied
   */
  virtual void upOpDim(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim);
};
}  // namespace finance
}  // namespace sgpp

#endif /* OPERATIONDELTALOGLINEARSTRETCHEDBOUNDARY_HPP */
