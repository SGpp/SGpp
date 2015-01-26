/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/algorithm/AlgorithmEvaluation.hpp>

#include <sgpp/base/basis/linear/noboundary/LinearBasis.hpp>

#include <sgpp/base/basis/linear/noboundary/operation/OperationEvalLinear.hpp>


namespace sg {
  namespace base {

    double OperationEvalLinear::eval(DataVector& alpha, std::vector<double>& point) {
      LinearBasis<unsigned int, unsigned int> base;
      AlgorithmEvaluation<LinearBasis<unsigned int, unsigned int> > AlgoEval(storage);

      return AlgoEval(base, point, alpha);
    }

  }
}

