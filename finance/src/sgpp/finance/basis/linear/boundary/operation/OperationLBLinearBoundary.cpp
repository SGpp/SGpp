/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include <sgpp/finance/basis/linear/boundary/operation/OperationLBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/DPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/DPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationLBLinearBoundary::OperationLBLinearBoundary(SGPP::base::GridStorage* storage) : SGPP::pde::StdUpDown(storage) {
    }

    OperationLBLinearBoundary::~OperationLBLinearBoundary() {
    }

    void OperationLBLinearBoundary::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // Dphi * phi
      DPhiPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<DPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLBLinearBoundary::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // Dphi * phi
      DPhiPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<DPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

  }
}
