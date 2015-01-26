// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>

#include <sgpp/base/basis/periodic/LinearPeriodicBasis.hpp>
#include <sgpp/base/basis/periodic/operation/OperationMultipleEvalPeriodic.hpp>



#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    void OperationMultipleEvalPeriodic::mult(DataVector& alpha, DataVector& result) {
      AlgorithmDGEMV<SLinearPeriodicBasis> op;
      LinearPeriodicBasis<unsigned int, unsigned int> base;

      op.mult(storage, base, alpha, this->dataset, result);

    }

    void OperationMultipleEvalPeriodic::multTranspose(DataVector& source, DataVector& result) {
      AlgorithmDGEMV<SLinearPeriodicBasis> op;
      LinearPeriodicBasis<unsigned int, unsigned int> base;

      op.mult_transposed(storage, base, source, this->dataset, result);
    }

  }
}