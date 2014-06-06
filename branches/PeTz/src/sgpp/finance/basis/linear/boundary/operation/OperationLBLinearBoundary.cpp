/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include "finance/basis/linear/boundary/operation/OperationLBLinearBoundary.hpp"

#include "finance/basis/linear/boundary/algorithm_sweep/DPhiPhiDownBBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/algorithm_sweep/DPhiPhiUpBBLinearBoundary.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg {
  namespace finance {

    OperationLBLinearBoundary::OperationLBLinearBoundary(sg::base::GridStorage* storage) : sg::pde::StdUpDown(storage) {
    }

    OperationLBLinearBoundary::~OperationLBLinearBoundary() {
    }

    void OperationLBLinearBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // Dphi * phi
      DPhiPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<DPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLBLinearBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // Dphi * phi
      DPhiPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<DPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

  }
}
