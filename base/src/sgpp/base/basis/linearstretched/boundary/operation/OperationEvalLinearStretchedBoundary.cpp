/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>


#include <sgpp/base/basis/linearstretched/boundary/operation/OperationEvalLinearStretchedBoundary.hpp>


namespace sg {
  namespace base {

    double OperationEvalLinearStretchedBoundary::eval(DataVector& alpha, std::vector<double>& point) {
      typedef std::vector<std::pair<size_t, double> > IndexValVector;

      IndexValVector vec;
      LinearStretchedBoundaryBasis<unsigned int, unsigned int> base;
      GetAffectedBasisFunctions<LinearStretchedBoundaryBasis<unsigned int, unsigned int> > ga(storage);

      ga(base, point, vec);

      double result = 0.0;

      for (IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++) {
        result += iter->second * alpha[iter->first];
      }

      return result;
    }

  }
}

