// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONMATRIXLTWODOTPOLYCLENSHAWCURTISBOUNDARY_HPP
#define OPERATIONMATRIXLTWODOTPOLYCLENSHAWCURTISBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace pde {

/**
 * Implements the standard L 2 scalar product on periodic grids
 *
 */
class OperationMatrixLTwoDotPolyClenshawCurtisBoundary : public sgpp::base::OperationMatrix {
 public:
  /**
   * Constructor
   *
   * @param grid  pointer to the sparse grid object
   */
  explicit OperationMatrixLTwoDotPolyClenshawCurtisBoundary(sgpp::base::Grid* grid);

  /**
   * Destructor
   */
  virtual ~OperationMatrixLTwoDotPolyClenshawCurtisBoundary();

  /**
  * Implementation of standard matrix multiplication
  *
  * @param alpha DataVector that is multiplied to the matrix
  * @param result DataVector into which the result of multiplication is stored
  */
  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

 protected:
  sgpp::base::Grid* grid;
  sgpp::base::ClenshawCurtisTable& clenshawCurtisTable;
};
}  // namespace pde
}  // namespace sgpp

#endif /* OPERATIONMATRIXLTWODOTPOLYCLENSHAWCURTISBOUNDARY_HPP */
