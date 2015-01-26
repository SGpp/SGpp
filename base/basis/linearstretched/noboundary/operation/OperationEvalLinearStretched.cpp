/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "base/algorithm/AlgorithmEvaluation.hpp"

#include "base/basis/linearstretched/noboundary/LinearStretchedBasis.hpp"

#include "base/basis/linearstretched/noboundary/operation/OperationEvalLinearStretched.hpp"


namespace sg {
  namespace base {

    double OperationEvalLinearStretched::eval(DataVector& alpha, std::vector<double>& point) {
      LinearStretchedBasis<unsigned int, unsigned int> base;
      AlgorithmEvaluation<LinearStretchedBasis<unsigned int, unsigned int> > AlgoEval(storage);

      return AlgoEval(base, point, alpha);
    }

  }
}

