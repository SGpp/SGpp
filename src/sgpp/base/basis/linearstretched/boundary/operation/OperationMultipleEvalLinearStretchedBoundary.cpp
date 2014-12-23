/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/algorithm/AlgorithmDGEMV.hpp"


#include "base/basis/linearstretched/boundary/operation/OperationMultipleEvalLinearStretchedBoundary.hpp"


namespace sg {
  namespace base {

    void OperationMultipleEvalLinearStretchedBoundary::mult(DataVector& alpha, DataVector& result) {
      AlgorithmDGEMV<SLinearStretchedBoundaryBase> op;
      LinearStretchedBoundaryBasis<unsigned int, unsigned int> base;

      op.mult(storage, base, alpha, this->dataset, result);
    }

    void OperationMultipleEvalLinearStretchedBoundary::multTranspose(DataVector& source, DataVector& result) {
      AlgorithmDGEMV<SLinearStretchedBoundaryBase> op;
      LinearStretchedBoundaryBasis<unsigned int, unsigned int> base;

      op.mult_transposed(storage, base, source, this->dataset, result);
    }

  }
}
