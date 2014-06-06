/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Sam Maurus (MA thesis)

#include "finance/basis/linear/noboundary/operation/pde/OperationHestonDLinear.hpp"

#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhidPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhidPhiUpBBLinear.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg {
  namespace finance {

    OperationHestonDLinear::OperationHestonDLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef) {
    }

    OperationHestonDLinear::~OperationHestonDLinear() {
    }

    void OperationHestonDLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonDLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonDLinear::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * dphi * dphi
      XdPhidPhiUpBBLinear func(this->storage);
      sg::base::sweep<XdPhidPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonDLinear::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * dphi * dphi
      XdPhidPhiDownBBLinear func(this->storage);
      sg::base::sweep<XdPhidPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

  }
}
