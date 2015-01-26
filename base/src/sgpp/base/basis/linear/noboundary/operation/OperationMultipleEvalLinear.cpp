/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp>
#include <sgpp/base/basis/linear/noboundary/LinearBasis.hpp>
#include <sgpp/base/basis/linear/noboundary/operation/OperationMultipleEvalLinear.hpp>


namespace sg {
  namespace base {

    void OperationMultipleEvalLinear::mult(DataVector& alpha, DataVector& result) {
      AlgorithmMultipleEvaluation<SLinearBase> op;
      LinearBasis<unsigned int, unsigned int> base;

      op.mult(storage, base, alpha, this->dataset, result);
    }

    void OperationMultipleEvalLinear::multTranspose(DataVector& alpha, DataVector& result) {
      AlgorithmMultipleEvaluation<SLinearBase> op;
      LinearBasis<unsigned int, unsigned int> base;

      op.mult_transpose(storage, base, alpha, this->dataset, result);
    }

  }
}
