// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/operation/hash/OperationEvalPrewavelet.hpp>

#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>



#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    float_t OperationEvalPrewavelet::eval(const DataVector& alpha,
                                          const DataVector& point) {
      typedef std::vector<std::pair<size_t, float_t> > IndexValVector;

      IndexValVector vec;
      PrewaveletBasis<unsigned int, unsigned int> base;
      GetAffectedBasisFunctions<PrewaveletBasis<unsigned int, unsigned int> > ga(
        storage);

      ga(base, point, vec);

      float_t result = 0.0;

      for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
        result += iter->second * alpha[iter->first];
      }

      return result;
    }

    float_t OperationEvalPrewavelet::test(const DataVector& alpha,
                                          const DataVector& data,
                                          const DataVector& classes) {
      return 0;
    }

    float_t OperationEvalPrewavelet::integrate(const DataVector& alpha) {
      float_t result = 0.0;

      for (size_t i = 0; i < storage->size(); i++) {
        float_t temp_result = 1;

        for (size_t d = 0; d < storage->dim(); d++) {
          GridStorage::index_type::level_type level;
          GridStorage::index_type::index_type index;
          (*storage)[i]->get(d, level, index);

          if (index != 1 && index != (unsigned int)((1 << level) - 1)) {
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

  }
}
