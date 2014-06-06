/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi(qic@in.tum.de)

#include "finance/basis/linear/noboundary/operation/pde/OperationLDLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XPhiPhiUpBBLinear.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg {
  namespace finance {

    OperationLDLinear::OperationLDLinear(sg::base::GridStorage* storage) : sg::pde::StdUpDown(storage) {
    }

    OperationLDLinear::~OperationLDLinear() {
    }

    void OperationLDLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // X * phi * phi
      XPhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<XPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationLDLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // X * phi * phi
      XPhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<XPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

  }
}
