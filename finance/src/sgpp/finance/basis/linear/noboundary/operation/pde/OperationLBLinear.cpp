/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include <sgpp/finance/basis/linear/noboundary/operation/pde/OperationLBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

namespace sg {
  namespace finance {

    OperationLBLinear::OperationLBLinear(sg::base::GridStorage* storage) : sg::pde::StdUpDown(storage) {
    }

    OperationLBLinear::~OperationLBLinear() {
    }

    void OperationLBLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // Dphi * phi
      DPhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<DPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationLBLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // Dphi * phi
      DPhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<DPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

  }
}
