/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#include <sgpp/base/basis/poly/operation/OperationEvalPoly.hpp>

#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    double OperationEvalPoly::eval(DataVector& alpha, std::vector<double>& point) {
      typedef std::vector<std::pair<size_t, double> > IndexValVector;

      IndexValVector vec;
      GetAffectedBasisFunctions<SPolyBase> ga(storage);

      ga(base, point, vec);

      double result = 0.0;

      for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
        result += iter->second * alpha[iter->first];
      }

      return result;
    }

  }
}
