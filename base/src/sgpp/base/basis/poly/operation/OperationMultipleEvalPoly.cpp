/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>

#include <sgpp/base/basis/poly/PolyBasis.hpp>
#include <sgpp/base/basis/poly/operation/OperationMultipleEvalPoly.hpp>



namespace sg {
  namespace base {

    void OperationMultipleEvalPoly::mult(DataVector& alpha, DataVector& result) {
      AlgorithmDGEMV<SPolyBase> op;

      op.mult(storage, base, alpha, this->dataset, result);
    }

    void OperationMultipleEvalPoly::multTranspose(DataVector& source, DataVector& result) {
      AlgorithmDGEMV<SPolyBase> op;

      op.mult_transposed(storage, base, source, this->dataset, result);
    }

  }
}
