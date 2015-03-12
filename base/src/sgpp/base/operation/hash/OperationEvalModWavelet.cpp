// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/operation/hash/OperationEvalModWavelet.hpp>

#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>



#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    float_t OperationEvalModWavelet::eval(DataVector& alpha, DataVector& point) {
      typedef std::vector<std::pair<size_t, float_t> > IndexValVector;

      IndexValVector vec;
      WaveletModifiedBasis<unsigned int, unsigned int> base;
      GetAffectedBasisFunctions<WaveletModifiedBasis<unsigned int, unsigned int> > ga(storage);

      ga(base, point, vec);

      float_t result = 0.0;

      for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
        result += iter->second * alpha[iter->first];
      }

      return result;
    }

  }
}
