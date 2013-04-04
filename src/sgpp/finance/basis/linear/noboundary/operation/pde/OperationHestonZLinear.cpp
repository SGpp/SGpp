/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Sam Maurus (MA thesis)

#include "finance/basis/linear/noboundary/operation/pde/OperationHestonZLinear.hpp"

#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg {
  namespace finance {

    OperationHestonZLinear::OperationHestonZLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef) {
    }

    OperationHestonZLinear::~OperationHestonZLinear() {
    }

    void OperationHestonZLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonZLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonZLinear::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * dphi * phi
      XdPhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<XdPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonZLinear::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * dphi * phi
      XdPhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<XdPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

  }
}

