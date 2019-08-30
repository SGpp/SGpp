// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/operation/hash/OperationEvalPrewavelet.hpp>

#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>



#include <sgpp/globaldef.hpp>

#include <utility>
#include <vector>


namespace sgpp {
namespace base {

double OperationEvalPrewavelet::eval(const DataVector& alpha,
                                      const DataVector& point) {
  typedef std::vector<std::pair<size_t, double> > IndexValVector;

  IndexValVector vec;
  PrewaveletBasis<unsigned int, unsigned int> base;
  GetAffectedBasisFunctions<PrewaveletBasis<unsigned int, unsigned int> > ga(storage);

  ga(base, point, vec);

  double result = 0.0;

  for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
    result += iter->second * alpha[iter->first];
  }

  return result;
}

double OperationEvalPrewavelet::test(const DataVector& alpha,
                                      const DataVector& data,
                                      const DataVector& classes) {
  return 0;
}

double OperationEvalPrewavelet::integrate(const DataVector& alpha) {
  double result = 0.0;

  for (size_t i = 0; i < storage.getSize(); i++) {
    double temp_result = 1;

    for (size_t d = 0; d < storage.getDimension(); d++) {
      level_t level;
      index_t index;
      storage[i].get(d, level, index);

      if (index != 1 && index != static_cast<unsigned int>((1 << level) - 1)) {
        temp_result = 0.0;
        break;
      } else if (level == 1) {
        temp_result *= 1.0 / 2.0;
      } else {
        temp_result *= 0.4 * (1.0 / (1 << level));
      }
    }

    result += alpha[i] * temp_result;
  }

  return result;
}

}  // namespace base
}  // namespace sgpp
