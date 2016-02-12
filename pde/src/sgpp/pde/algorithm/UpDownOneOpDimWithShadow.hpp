// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef UPDOWNONEOPDIMWITHSHADOW_HPP
#define UPDOWNONEOPDIMWITHSHADOW_HPP

#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/operation/hash/OperationMatrix.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace pde {

/**
 * Implements the Up/Down scheme with one dimension with a special operation. Before the actual
 * operation starts, all grid points from the shadow storage are copied into the actual grid.
 * After the calculation, all shadow grid points are deleted, the result vector is adapted
 * accordingly.
 *
 */
class UpDownOneOpDimWithShadow : public SGPP::base::OperationMatrix {
 public:
  /**
   * Constructor
   *
   * @param storage the grid's SGPP::base::GridStorage object
   * @param shadowStorage shadow storage
   */
  UpDownOneOpDimWithShadow(SGPP::base::GridStorage* storage,
                           SGPP::base::GridStorage* shadowStorage);

  /**
   * Destructor
   */
  virtual ~UpDownOneOpDimWithShadow();

  virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

 protected:
  typedef SGPP::base::GridStorage::grid_iterator grid_iterator;

  /// Pointer to the grid's storage object
  SGPP::base::GridStorage* storage;
  SGPP::base::GridStorage* shadowStorage;

  /**
   * Recursive procedure for updown().
   *
   * @param dim the current dimension
   * @param op_dim the dimension in which a special operation is applied
   * @param alpha vector of coefficients
   * @param result vector to store the results in
   */
  void updown(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim,
              size_t op_dim);

  /**
   * This functions adds all grid points of the shadow storage into the actual grid.
   */
  void expandGrid();

  /**
   * Removes the previously added shadow grid points from the actual grid. Thus, the
   * grid is the same shape as before the call of mult.
   */
  void shrinkGrid();

  /**
   * All calculations for op_dim.
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

#endif /* UPDOWNONEOPDIMWITHSHADOW_HPP */
