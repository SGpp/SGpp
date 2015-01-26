/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include <sgpp/finance/basis/linear/boundary/operation/OperationLELinearBoundary.hpp>


#include <sgpp/pde/basis/linear/boundary/operation/DowndPhidPhiBBIterativeLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/operation/UpdPhidPhiBBIterativeLinearBoundary.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationLELinearBoundary::OperationLELinearBoundary(SGPP::base::GridStorage* storage) : SGPP::pde::StdUpDown(storage) {
    }

    OperationLELinearBoundary::~OperationLELinearBoundary() {
    }

    void OperationLELinearBoundary::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // Dphi * dphi
      SGPP::pde::UpdPhidPhiBBIterativeLinearBoundary myUp(this->storage);
      myUp(alpha, result, dim);
    }

    void OperationLELinearBoundary::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // Dphi * dphi
      SGPP::pde::DowndPhidPhiBBIterativeLinearBoundary myDown(this->storage);
      myDown(alpha, result, dim);
    }

  }
}
