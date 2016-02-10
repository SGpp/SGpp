// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONLTWODOTPRODUCTVECTORIZEDLINEARBOUNDARYOCL_HPP
#define OPERATIONLTWODOTPRODUCTVECTORIZEDLINEARBOUNDARYOCL_HPP

#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/parallel/pde/basis/common/OCLPDEKernels.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace parallel {

/**
 * Implements the standard L 2 scalar product on linear boundary grids
 *
 */
class OperationLTwoDotProductVectorizedLinearBoundaryOCL: public
  SGPP::base::OperationMatrix {

 private:
  SGPP::base::GridStorage* storage;
  SGPP::base::DataMatrix* level_;
  SGPP::base::DataMatrix* level_int_;
  SGPP::base::DataMatrix* index_;
  double* lcl_q;
  OCLPDEKernels OCLPDEKernelsHandle;

  void mult_dirichlet(SGPP::base::DataVector& alpha,
                      SGPP::base::DataVector& result);

 public:
  /**
   * Constructor
   *
   * @param storage the grid's SGPP::base::GridStorage object
   */
  OperationLTwoDotProductVectorizedLinearBoundaryOCL(SGPP::base::GridStorage*
      storage);

  /**
   * Destructor
   */
  virtual ~OperationLTwoDotProductVectorizedLinearBoundaryOCL();

 protected:
  virtual void mult(SGPP::base::DataVector& alpha,
                    SGPP::base::DataVector& result);
};

}
}

#endif /* OPERATIONLTWODOTPRODUCTVECTORIZEDLINEARBOUNDARYOCL_HPP */