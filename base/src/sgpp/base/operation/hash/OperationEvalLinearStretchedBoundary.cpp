// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>


#include <sgpp/base/operation/hash/OperationEvalLinearStretchedBoundary.hpp>


#include <sgpp/globaldef.hpp>

#include <utility>
#include <vector>


namespace SGPP {
namespace base {

float_t OperationEvalLinearStretchedBoundary::eval(const DataVector& alpha,
    const DataVector& point) {
  typedef std::vector<std::pair<size_t, float_t> > IndexValVector;

  IndexValVector vec;
  LinearStretchedBoundaryBasis<unsigned int, unsigned int> base;
  GetAffectedBasisFunctions <
  LinearStretchedBoundaryBasis<unsigned int, unsigned int> >
  ga(storage);

  ga(base, point, vec);

  float_t result = 0.0;

  for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
    result += iter->second * alpha[iter->first];
  }

  return result;
}

}  // namespace base
}  // namespace SGPP
