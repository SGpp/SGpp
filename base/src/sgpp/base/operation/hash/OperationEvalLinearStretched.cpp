// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmEvaluation.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearStretchedBasis.hpp>

#include <sgpp/base/operation/hash/OperationEvalLinearStretched.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    double OperationEvalLinearStretched::eval(DataVector& alpha, std::vector<double>& point) {
      LinearStretchedBasis<unsigned int, unsigned int> base;
      AlgorithmEvaluation<LinearStretchedBasis<unsigned int, unsigned int> > AlgoEval(storage);

      return AlgoEval(base, point, alpha);
    }

  }
}
