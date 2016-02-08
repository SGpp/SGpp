// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLELINEARBOUNDARY_HPP
#define OPERATIONLELINEARBOUNDARY_HPP

#include <sgpp/pde/algorithm/StdUpDown.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace finance {

/**
 * Implements the \f$(d\phi_i(x),d\phi_j(x))\f$ operator on linear boundary grids
 *
 */
class OperationLELinearBoundary: public SGPP::pde::StdUpDown {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's SGPP::base::GridStorage object
   */
  OperationLELinearBoundary(SGPP::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~OperationLELinearBoundary();

 protected:
  /**
   * Up-step in dimension <i>dim</i> for \f$(d\phi_i(x),d\phi_j(x))\f$.
   * Applies the up-part of the one-dimensional mass matrix in one dimension.
   * Computes \f[\int_{x=0}^1  d\phi_i(x) d\phi_j(x) dx.\f]
   *
   * @param dim dimension in which to apply the up-part
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  virtual void up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                  size_t dim);

  /**
   * Down-step in dimension <i>dim</i> for \f$(d\phi_i(x),d\phi_j(x))\f$.
   * Applies the down-part of the one-dimensional mass matrix in one dimension.
   * Computes \f[\int_{x=0}^1  d\phi_i(x) d\phi_j(x) dx.\f]
   *
   * @param dim dimension in which to apply the down-part
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  virtual void down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                    size_t dim);
};

}
}

#endif /* OPERATIONLELINEARBOUNDARY_HPP */