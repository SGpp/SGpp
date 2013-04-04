/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/algorithm/GetAffectedBasisFunctions.hpp"


#include "base/basis/linear/boundary/operation/OperationEvalLinearBoundary.hpp"


namespace sg {
  namespace base {

    double OperationEvalLinearBoundary::eval(DataVector& alpha, std::vector<double>& point) {
      typedef std::vector<std::pair<size_t, double> > IndexValVector;

      IndexValVector vec;
      LinearBoundaryBasis<unsigned int, unsigned int> base;
      GetAffectedBasisFunctions<LinearBoundaryBasis<unsigned int, unsigned int> > ga(storage);

      ga(base, point, vec);

      double result = 0.0;

      for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
        result += iter->second * alpha[iter->first];
      }

      return result;
    }

  }
}

