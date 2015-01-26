/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/basis/modpoly/ModifiedPolyBasis.hpp>
#include <sgpp/base/basis/modwavelet/operation/OperationMultipleEvalModWavelet.hpp>

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>



namespace sg {
  namespace base {

    void OperationMultipleEvalModWavelet::mult(DataVector& alpha, DataVector& result) {
      AlgorithmDGEMV<SModWaveletBase> op;
      ModifiedWaveletBasis<unsigned int, unsigned int> base;

      op.mult(storage, base, alpha, this->dataset, result);
    }

    void OperationMultipleEvalModWavelet::multTranspose(DataVector& source, DataVector& result) {
      AlgorithmDGEMV<SModWaveletBase> op;
      ModifiedWaveletBasis<unsigned int, unsigned int> base;

      op.mult_transposed(storage, base, source, this->dataset, result);
    }

  }
}
