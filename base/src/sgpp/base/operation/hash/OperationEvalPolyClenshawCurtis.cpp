// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationEvalPolyClenshawCurtis.hpp>

#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <utility>
#include <vector>

namespace sgpp {
namespace base {

double OperationEvalPolyClenshawCurtis::eval(const DataVector& alpha, const DataVector& point) {
  typedef std::vector<std::pair<size_t, double> > IndexValVector;

  IndexValVector vec;
  GetAffectedBasisFunctions<SPolyClenshawCurtisBase> ga(storage);

  ga(base, point, vec);

  double result = 0.0;

  for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
    result += iter->second * alpha[iter->first];
  }

  return result;
}

}  // namespace base
}  // namespace sgpp
