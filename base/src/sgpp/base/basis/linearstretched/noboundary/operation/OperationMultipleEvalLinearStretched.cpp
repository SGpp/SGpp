/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/algorithm/AlgorithmDGEMV.hpp"
#include "base/basis/linearstretched/noboundary/LinearStretchedBasis.hpp"
#include "base/basis/linearstretched/noboundary/operation/OperationMultipleEvalLinearStretched.hpp"


namespace sg {
  namespace base {

    void OperationMultipleEvalLinearStretched::mult(DataVector& alpha, DataVector& result) {
      AlgorithmDGEMV<SLinearStretchedBase> op;
      LinearStretchedBasis<unsigned int, unsigned int> base;

      op.mult(storage, base, alpha, this->dataset, result);
    }

    void OperationMultipleEvalLinearStretched::multTranspose(DataVector& source, DataVector& result) {
      AlgorithmDGEMV<SLinearStretchedBase> op;
      LinearStretchedBasis<unsigned int, unsigned int> base;

      op.mult_transposed(storage, base, source, this->dataset, result);
    }

  }
}
