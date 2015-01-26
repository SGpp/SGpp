/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>

#include <sgpp/base/basis/modlinear/ModifiedLinearBasis.hpp>
#include <sgpp/base/basis/modlinear/operation/OperationMultipleEvalModLinear.hpp>



#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    void OperationMultipleEvalModLinear::mult(DataVector& alpha, DataVector& result) {
      AlgorithmDGEMV<SModLinearBase> op;
      ModifiedLinearBasis<unsigned int, unsigned int> base;

      op.mult(storage, base, alpha, this->dataset, result);
    }

    void OperationMultipleEvalModLinear::multTranspose(DataVector& source, DataVector& result) {
      AlgorithmDGEMV<SModLinearBase> op;
      ModifiedLinearBasis<unsigned int, unsigned int> base;

      op.mult_transposed(storage, base, source, this->dataset, result);
    }

  }
}
