/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include <sgpp/finance/basis/linear/noboundary/operation/pde/OperationLELinear.hpp>

//#include "basis/linear/noboundary/algorithm_sweep/DPhidPhiDownBBLinear.hpp"
//#include "basis/linear/noboundary/algorithm_sweep/DPhidPhiUpBBLinear.hpp"
#include <sgpp/pde/basis/linear/noboundary/DowndPhidPhiBBIterativeLinear.hpp>
#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationLELinear::OperationLELinear(SGPP::base::GridStorage* storage) : SGPP::pde::StdUpDown(storage) {
    }

    OperationLELinear::~OperationLELinear() {
    }

    void OperationLELinear::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {

    }

    void OperationLELinear::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // Dphi * dphi
      SGPP::pde::DowndPhidPhiBBIterativeLinear myDown(this->storage);
      myDown(alpha, result, dim);
    }

  }
}
