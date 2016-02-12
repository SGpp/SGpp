// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef UPDOWNONEOPDIM_HPP
#define UPDOWNONEOPDIM_HPP

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#ifndef TASKS_PARALLEL_UPDOWN
#define TASKS_PARALLEL_UPDOWN 4
#endif

#include <sgpp/globaldef.hpp>

#include <vector>

namespace SGPP {
namespace pde {

/**
 * Implements the Up/Down scheme with one dimension with a special operation
 *
 * Parallelization with OpenMP 2 / 3 is supported!
 *
 */
class UpDownOneOpDim : public SGPP::base::OperationMatrix {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's SGPP::base::GridStorage object
   * @param coef reference to a SGPP::base::DataVector object that contains the bilinear form's
   * constant coefficients; one per dimension
   */
  UpDownOneOpDim(SGPP::base::GridStorage* storage, SGPP::base::DataVector& coef);

  /**
   * Constructor
   *
   * @param storage the grid's SGPP::base::GridStorage object
   */
  explicit UpDownOneOpDim(SGPP::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~UpDownOneOpDim();

  virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

  /**
   * This functions provides the same functionality as the normal mult routine.
   * However, it doesn't set up an OpenMP task initialization as the mult routine.
   * This method has to be called within a OpenMP task parallelized region.
   *
   * Using this function is useful in following case: Assuming the solver of a certain
   * requires several operators in the space discretization (e.g. Black Scholes Equations)
   * this method can be used to parallelize their calculation which might results results
   * in a better parallel efficiency on systems with 4 or more cores hence fewer barriers
   * are needed.
   *
   * For full calculation in of mult serval number of up/downs are needed. This number
   * is equal to the number of the grid's dimensions. All different steps can be executed
   * in parallel, so here only one up/Down is executed, identified by its special dimension.
   *
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   * @param operationDim Dimension in which the special operator is applied
   */
  void multParallelBuildingBlock(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                                 size_t operationDim);

 protected:
  typedef SGPP::base::GridStorage::grid_iterator grid_iterator;

  /// Pointer to the grid's storage object
  SGPP::base::GridStorage* storage;
  /// Pointer to the SGPP::base::DataVector of the coefs
  SGPP::base::DataVector* coefs;
  /// algorithmic dimensions, operator is applied in this dimensions
  const std::vector<size_t> algoDims;
  /// number of algorithmic dimensions
  const size_t numAlgoDims_;
  /// max number of parallel stages (dimension recursive calls)
  static const size_t maxParallelDims_ = TASKS_PARALLEL_UPDOWN;

  /**
   * Recursive procedure for updown(), parallel version using OpenMP 3
   *
   * @param dim the current dimension
   * @param op_dim the dimension in which a special operation is applied
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  void updown(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim,
              size_t op_dim);

  /**
   * All calculations for gradient_dim, parallel version using OpenMP 3
   *
   * @param alpha the coefficients of the grid points
   * @param result the result of the operations
   * @param dim the current dimension in the recursion
   * @param op_dim the dimension in that a special operation is applied
   */
  virtual void specialOP(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim,
                         size_t op_dim);

  /**
   * std 1D up operation
   *
   * @param dim dimension in which to apply the up-part
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  virtual void up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) = 0;

  /**
   * std 1D down operation
   *
   * @param dim dimension in which to apply the down-part
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  virtual void down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) = 0;

  /**
   * special 1D down operation that is only executed in one direction
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that down-Gradient is applied
   */
  virtual void downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                         size_t dim) = 0;

  /**
   * special 1D up operation that is only executed in one direction
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that up-Gradient is applied
   */
  virtual void upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                       size_t dim) = 0;
};
}  // namespace pde
}  // namespace SGPP

#endif /* UPDOWNONEOPDIM_HPP */
