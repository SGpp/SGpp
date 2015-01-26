/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Alexander Heinecke (alexander.heinecke@mytum.de)

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>

#include <sgpp/base/basis/modbspline/operation/OperationMultipleEvalModBspline.hpp>



namespace sg {
  namespace base {

    void OperationMultipleEvalModBspline::mult(DataVector& alpha, DataVector& result) {
      AlgorithmDGEMV<SModBsplineBase> op;

      op.mult(storage, base, alpha, this->dataset, result);
    }

    void OperationMultipleEvalModBspline::multTranspose(DataVector& source, DataVector& result) {
      AlgorithmDGEMV<SModBsplineBase> op;

      op.mult_transposed(storage, base, source, this->dataset, result);
    }

  }
}
