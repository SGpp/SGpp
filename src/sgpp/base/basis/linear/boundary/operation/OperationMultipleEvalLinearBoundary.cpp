/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/algorithm/AlgorithmDGEMV.hpp"


#include "base/basis/linear/boundary/operation/OperationMultipleEvalLinearBoundary.hpp"


namespace sg {
  namespace base {

    void OperationMultipleEvalLinearBoundary::mult(DataVector& alpha, DataVector& result) {
      AlgorithmDGEMV<SLinearBoundaryBase> op;
      LinearBoundaryBasis<unsigned int, unsigned int> base;

      op.mult(storage, base, alpha, this->dataset, result);
    }

    void OperationMultipleEvalLinearBoundary::multTranspose(DataVector& source, DataVector& result) {
      AlgorithmDGEMV<SLinearBoundaryBase> op;
      LinearBoundaryBasis<unsigned int, unsigned int> base;

      op.mult_transposed(storage, base, source, this->dataset, result);
    }

  }
}
