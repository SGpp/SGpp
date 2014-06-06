/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include "finance/basis/linear/noboundary/operation/pde/OperationLELinear.hpp"

//#include "basis/linear/noboundary/algorithm_sweep/DPhidPhiDownBBLinear.hpp"
//#include "basis/linear/noboundary/algorithm_sweep/DPhidPhiUpBBLinear.hpp"
#include "pde/basis/linear/noboundary/DowndPhidPhiBBIterativeLinear.hpp"
#include "base/algorithm/sweep.hpp"

namespace sg {
  namespace finance {

    OperationLELinear::OperationLELinear(sg::base::GridStorage* storage) : sg::pde::StdUpDown(storage) {
    }

    OperationLELinear::~OperationLELinear() {
    }

    void OperationLELinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {

    }

    void OperationLELinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // Dphi * dphi
      sg::pde::DowndPhidPhiBBIterativeLinear myDown(this->storage);
      myDown(alpha, result, dim);
    }

  }
}
