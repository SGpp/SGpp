/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include <sgpp/finance/basis/linear/boundary/operation/OperationLDLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/XPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

namespace sg {
  namespace finance {

    OperationLDLinearBoundary::OperationLDLinearBoundary(sg::base::GridStorage* storage) : sg::pde::StdUpDown(storage) {
    }

    OperationLDLinearBoundary::~OperationLDLinearBoundary() {
    }

    void OperationLDLinearBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // X * phi * phi
      XPhiPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<XPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLDLinearBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // X * phi * phi
      XPhiPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<XPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

  }
}
