// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLTWODOTPRODUCTLINEARSTRETCHED_HPP
#define OPERATIONLTWODOTPRODUCTLINEARSTRETCHED_HPP

#include <sgpp/pde/algorithm/StdUpDown.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Implements the standard L 2 scalar product on linear grids (no boundaries)
 *
 */
class OperationLTwoDotProductLinearStretched : public StdUpDown {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's sgpp::base::GridStorage object
   */
  explicit OperationLTwoDotProductLinearStretched(sgpp::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~OperationLTwoDotProductLinearStretched();

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
};
}  // namespace pde
}  // namespace sgpp

#endif /* OPERATIONLTWODOTPRODUCTLINEARSTRETCHED_HPP */
