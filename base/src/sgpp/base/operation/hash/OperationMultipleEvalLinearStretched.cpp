// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>
#include <sgpp/base/basis/linearstretched/noboundary/LinearStretchedBasis.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinearStretched.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
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