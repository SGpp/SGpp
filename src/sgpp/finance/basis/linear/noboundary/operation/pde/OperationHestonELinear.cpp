/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Sam Maurus (MA thesis)

#include "finance/basis/linear/noboundary/operation/pde/OperationHestonELinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp"

#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg {
  namespace finance {

    OperationHestonELinear::OperationHestonELinear(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef) {
    }

    OperationHestonELinear::~OperationHestonELinear() {
    }

    void OperationHestonELinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);

    }

    void OperationHestonELinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);

    }

    void OperationHestonELinear::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<DPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonELinear::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<DPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

  }
}
