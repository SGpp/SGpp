// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef UPDOWNTWOOPDIMS_HPP
#define UPDOWNTWOOPDIMS_HPP

#include <vector>

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#ifndef TASKS_PARALLEL_UPDOWN
#define TASKS_PARALLEL_UPDOWN 4
#endif

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace pde {

/**
 * Implements the Up/Down scheme with two dimensions with special operations: i,j
 *
 * Parallelization with OpenMP 2 / 3 is supported!
 *
 * Only symmetric operations are support --> only
 * the "left lower triangular matrix", i <= j, is calculated, please
 * keep that in mind when designing the coefficient vector:
 * the non-diagonal elements must be multiplied by 2
 * before executing this Up/down scheme!
 *
 */
class UpDownTwoOpDims: public SGPP::base::OperationMatrix {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's SGPP::base::GridStorage object
   * @param coef vector that contains the constant coefficients of this operation
   */
  UpDownTwoOpDims(SGPP::base::GridStorage* storage, SGPP::base::DataMatrix& coef);

  /**
   * Constructor
   *
   * @param storage the grid's SGPP::base::GridStorage object
   */
  UpDownTwoOpDims(SGPP::base::GridStorage* storage);

  /**
   * Destructor
   */
  virtual ~UpDownTwoOpDims();

  virtual void mult(SGPP::base::DataVector& alpha,
                    SGPP::base::DataVector& result);

  /**
   * this functions provides the same functionality as the normal mult routine.
   * However, it doesn't set up an OpenMP task initialization as the mult routine.
   * This method has to be called within a OpenMP task parallelized region.
   *
   * Using this function is useful in the following case: Assuming the solver of a certain Equation
   * requires several operators in the space discretization (e.g. Black Scholes Equations)
   * this method can be used to parallelize their calculation which might result
   * in a better parallel efficiency on systems with 4 or more cores hence fewer barriers
   * are needed.
   *
   * For a full calculation of this operator, in mult serval number of up/downs are needed. This number
   * is equal to the square of the number of the grid's dimensions. All different steps can be executed
   * in parallel. Here only one up/Down is executed, identified by its special dimensions.
   *
   * Attention: A symmetric operator is assumed: This method only start a calculation if
   * operationDimTwo is less or equal operationDimOne
   *
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   * @param operationDimOne Dimension in which the first special operator is applied
   * @param operationDimTwo Dimension in which the second special operator is applied
   */
  void multParallelBuildingBlock(SGPP::base::DataVector& alpha,
                                 SGPP::base::DataVector& result, size_t operationDimOne, size_t operationDimTwo);

 protected:
  typedef SGPP::base::GridStorage::grid_iterator grid_iterator;

  /// Pointer to the grid's storage object
  SGPP::base::GridStorage* storage;
  /// Pointer to the coefficients of this bilinear form
  SGPP::base::DataMatrix* coefs;
  /// algorithmic dimensions, operator is applied in this dimensions
  const std::vector<size_t> algoDims;
  /// number of algorithmic dimensions
  const size_t numAlgoDims_;
  /// max number of parallel stages (dimension recursive calls)
  static const size_t maxParallelDims_ = TASKS_PARALLEL_UPDOWN;

  /**
   * Recursive procedure for updown, parallel version using OpenMP 3
   *
   * @param dim the current dimension
   * @param op_dim_one the dimension in which to use the first gradient
   * @param op_dim_two the dimension in which to use the second gradient
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  void updown(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
              size_t dim, size_t op_dim_one, size_t op_dim_two);

  /**
   * All calculations for gradient, parallel version using OpenMP 3
   *
   * @param alpha the coefficients of the grid points
   * @param result the result of the operations
   * @param dim the current dimension in the recursion
   * @param op_dim_one the dimension in which to use the first gradient
   * @param op_dim_two the dimension in which to use the second gradient
   */
  void specialOpOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                    size_t dim, size_t op_dim_one, size_t op_dim_two);

  /**
   * All calculations for gradient, Part 2, parallel version using OpenMP 3
   *
   * @param alpha the coefficients of the grid points
   * @param result the result of the operations
   * @param dim the current dimension in the recursion
   * @param op_dim_one the dimension in which to use the first gradient
   * @param op_dim_two the dimension in which to use the second gradient
   */
  void specialOpTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                    size_t dim, size_t op_dim_one, size_t op_dim_two);

  /**
   * if the current dimension is equal to the both special operation dimensions
   *
   * @param alpha the coefficients of the grid points
   * @param result the result of the operations
   * @param dim the current dimension in the recursion
   * @param op_dim_one the dimension in which to use the first gradient
   * @param op_dim_two the dimension in which to use the second gradient
   */
  void specialOpOneAndOpTwo(SGPP::base::DataVector& alpha,
                            SGPP::base::DataVector& result, size_t dim, size_t op_dim_one,
                            size_t op_dim_two);

  /**
   * Up-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
   * Applies the up-part of the one-dimensional mass matrix in one dimension.
   * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i < l_j} \alpha_j \phi_j(x) dx.\f]
   *
   * @param dim dimension in which to apply the up-part
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  virtual void up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                  size_t dim) = 0;

  /**
   * Down-step in dimension <i>dim</i> for \f$(\phi_i(x),\phi_j(x))_{L_2}\f$.
   * Applies the down-part of the one-dimensional mass matrix in one dimension.
   * Computes \f[\int_{x=0}^1  \phi_i(x) \sum_{j, l_i\geq l_j} \alpha_j \phi_j(x) dx.\f]
   *
   * @param dim dimension in which to apply the down-part
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  virtual void down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                    size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to i
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that down-Gradient is applied
   */
  virtual void downOpDimOne(SGPP::base::DataVector& alpha,
                            SGPP::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to i
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that up-Gradient is applied
   */
  virtual void upOpDimOne(SGPP::base::DataVector& alpha,
                          SGPP::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D down if the current dim is equal to j
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that down-Gradient is applied
   */
  virtual void downOpDimTwo(SGPP::base::DataVector& alpha,
                            SGPP::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D up if the current dim is equal to j
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that up-Gradient is applied
   */
  virtual void upOpDimTwo(SGPP::base::DataVector& alpha,
                          SGPP::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D down, if the current dim is equal to i and j
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that down-Gradient is applied
   */
  virtual void downOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha,
                                       SGPP::base::DataVector& result, size_t dim) = 0;

  /**
   * 1D up, if the current dim is equal to i and j
   *
   * @param alpha the coefficients of the gridpoints
   * @param result vector with the result of this operation
   * @param dim the dimension in that up-Gradient is applied
   */
  virtual void upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha,
                                     SGPP::base::DataVector& result, size_t dim) = 0;
};

}
}

#endif /* UPDOWNTWOOPDIMS_HPP */