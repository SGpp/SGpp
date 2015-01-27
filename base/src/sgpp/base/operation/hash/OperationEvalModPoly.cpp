// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationEvalModPoly.hpp>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    double OperationEvalModPoly::eval(DataVector& alpha, std::vector<double>& point) {
      typedef std::vector<std::pair<size_t, double> > IndexValVector;

      IndexValVector vec;
      GetAffectedBasisFunctions<SPolyModifiedBase> ga(storage);

      ga(base, point, vec);

      double result = 0.0;

      for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
        result += iter->second * alpha[iter->first];
      }

      return result;
    }

  }
}