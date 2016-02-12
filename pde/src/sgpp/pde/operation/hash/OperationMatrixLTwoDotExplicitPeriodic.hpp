// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMATRIXLTWODOTEXPLICITPERIODIC_HPP
#define OPERATIONMATRIXLTWODOTEXPLICITPERIODIC_HPP

#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace pde {

/**
 * Explicit representation of the matrix \f$(\Phi_i,\Phi_j)_{L2}\f$ for a sparse grid
 */
class OperationMatrixLTwoDotExplicitPeriodic : public SGPP::base::OperationMatrix {
 public:
  /**
   * Constructor that uses a external matrix pointer to construct the matrix,
   *
   * @param m pointer to datamatrix of size (number of grid point) x (number of grid points)
   * @param grid the sparse grid
   */
  OperationMatrixLTwoDotExplicitPeriodic(SGPP::base::DataMatrix* m, SGPP::base::Grid* grid);
  /**
   * Constructor that creates an own matrix
   *
   * @param grid the sparse grid
   */
  explicit OperationMatrixLTwoDotExplicitPeriodic(SGPP::base::Grid* grid);

  /**
   * Destructor
   */
  virtual ~OperationMatrixLTwoDotExplicitPeriodic();

  /**
   * Implementation of standard matrix multiplication
   *
   * @param alpha DataVector that is multiplied to the matrix
   * @param result DataVector into which the result of multiplication is stored
   */
  virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

 private:
  /**
   * This method is used by both constructors to build the matrix
   */
  void buildMatrix(SGPP::base::Grid* grid);

  SGPP::base::DataMatrix* m_;
  bool ownsMatrix_;
};

}  // namespace pde
}  // namespace SGPP

#endif /* OPERATIONMATRIXLTWODOTEXPLICITPERIODIC_HPP */
