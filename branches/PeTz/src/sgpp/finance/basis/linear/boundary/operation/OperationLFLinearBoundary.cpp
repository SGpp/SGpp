/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include "finance/basis/linear/boundary/operation/OperationLFLinearBoundary.hpp"

#include "finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiDownBBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/algorithm_sweep/XdPhiPhiUpBBLinearBoundary.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg {
  namespace finance {

    OperationLFLinearBoundary::OperationLFLinearBoundary(sg::base::GridStorage* storage) : sg::pde::StdUpDown(storage) {
    }

    OperationLFLinearBoundary::~OperationLFLinearBoundary() {
    }

    void OperationLFLinearBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // X * dphi * phi
      XdPhiPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<XdPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLFLinearBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // X * dphi * phi
      XdPhiPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<XdPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

  }
}
