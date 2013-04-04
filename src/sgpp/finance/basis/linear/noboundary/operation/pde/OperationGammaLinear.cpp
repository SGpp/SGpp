/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "finance/basis/linear/noboundary/operation/pde/OperationGammaLinear.hpp"

#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XPhidPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XPhidPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/XdPhiPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/SqXdPhidPhiUpBBLinear.hpp"

#include "base/algorithm/sweep.hpp"

#include <iostream>

namespace sg {
  namespace finance {

    OperationGammaLinear::OperationGammaLinear(sg::base::GridStorage* storage, sg::base::DataMatrix& coef) : sg::pde::UpDownTwoOpDims(storage, coef) {
    }

    OperationGammaLinear::~OperationGammaLinear() {
    }

    void OperationGammaLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinear::upOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * phi * dphi
      XPhidPhiUpBBLinear func(this->storage);
      sg::base::sweep<XPhidPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinear::downOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * phi * dphi
      XPhidPhiDownBBLinear func(this->storage);
      sg::base::sweep<XPhidPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinear::upOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * dphi * phi
      XdPhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<XdPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinear::downOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * dphi * phi
      XdPhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<XdPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinear::upOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x^2 * dphi * dphi
      SqXdPhidPhiUpBBLinear func(this->storage);
      sg::base::sweep<SqXdPhidPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinear::downOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x^2 * dphi * dphi
      SqXdPhidPhiDownBBLinear func(this->storage);
      sg::base::sweep<SqXdPhidPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

  }
}
