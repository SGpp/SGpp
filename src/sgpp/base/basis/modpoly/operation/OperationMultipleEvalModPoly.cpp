/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/algorithm/AlgorithmDGEMV.hpp"

#include "base/basis/modpoly/ModifiedPolyBasis.hpp"
#include "base/basis/modpoly/operation/OperationMultipleEvalModPoly.hpp"



namespace sg {
  namespace base {

    void OperationMultipleEvalModPoly::mult(DataVector& alpha, DataVector& result) {
      AlgorithmDGEMV<SModPolyBase> op;

      op.mult(storage, base, alpha, this->dataset, result);
    }

    void OperationMultipleEvalModPoly::multTranspose(DataVector& source, DataVector& result) {
      AlgorithmDGEMV<SModPolyBase> op;

      op.mult_transposed(storage, base, source, this->dataset, result);
    }

  }
}
