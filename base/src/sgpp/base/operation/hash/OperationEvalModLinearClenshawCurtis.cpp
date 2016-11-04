// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/AlgorithmEvaluation.hpp>

#include <sgpp/base/operation/hash/common/basis/LinearModifiedClenshawCurtisBasis.hpp>

#include <sgpp/base/operation/hash/OperationEvalModLinearClenshawCurtis.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

double OperationEvalModLinearClenshawCurtis::eval(const DataVector& alpha,
                                                  const DataVector& point) {
  LinearModifiedClenshawCurtisBasis<unsigned int, unsigned int> base;
  AlgorithmEvaluation<LinearModifiedClenshawCurtisBasis<unsigned int, unsigned int> > AlgoEval(
      storage);

  return AlgoEval(base, point, alpha);
}

}  // namespace base
}  // namespace sgpp
