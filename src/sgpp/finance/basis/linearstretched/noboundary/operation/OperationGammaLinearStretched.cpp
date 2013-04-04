/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "finance/basis/linearstretched/noboundary/operation/OperationGammaLinearStretched.hpp"

#include "pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp"
#include "pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp"

#include "finance/basis/linearstretched/noboundary/algorithm_sweep/XPhidPhiDownBBLinearStretched.hpp"
#include "finance/basis/linearstretched/noboundary/algorithm_sweep/XPhidPhiUpBBLinearStretched.hpp"

#include "finance/basis/linearstretched/noboundary/algorithm_sweep/XdPhiPhiDownBBLinearStretched.hpp"
#include "finance/basis/linearstretched/noboundary/algorithm_sweep/XdPhiPhiUpBBLinearStretched.hpp"

#include "finance/basis/linearstretched/noboundary/algorithm_sweep/SqXdPhidPhiDownBBLinearStretched.hpp"
#include "finance/basis/linearstretched/noboundary/algorithm_sweep/SqXdPhidPhiUpBBLinearStretched.hpp"

#include "base/algorithm/sweep.hpp"

#include <iostream>

namespace sg {
  namespace finance {

    OperationGammaLinearStretched::OperationGammaLinearStretched(sg::base::GridStorage* storage, sg::base::DataMatrix& coef) : sg::pde::UpDownTwoOpDims(storage, coef) {
    }

    OperationGammaLinearStretched::~OperationGammaLinearStretched() {
    }

    void OperationGammaLinearStretched::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiUpBBLinearStretched func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiUpBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinearStretched::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiDownBBLinearStretched func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiDownBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinearStretched::upOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * phi * dphi
      XPhidPhiUpBBLinearStretched func(this->storage);
      sg::base::sweep<XPhidPhiUpBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinearStretched::downOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * phi * dphi
      XPhidPhiDownBBLinearStretched func(this->storage);
      sg::base::sweep<XPhidPhiDownBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinearStretched::upOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * dphi * phi
      XdPhiPhiUpBBLinearStretched func(this->storage);
      sg::base::sweep<XdPhiPhiUpBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinearStretched::downOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x * dphi * phi
      XdPhiPhiDownBBLinearStretched func(this->storage);
      sg::base::sweep<XdPhiPhiDownBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinearStretched::upOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x^2 * dphi * dphi
      SqXdPhidPhiUpBBLinearStretched func(this->storage);
      sg::base::sweep<SqXdPhidPhiUpBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationGammaLinearStretched::downOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // x^2 * dphi * dphi
      SqXdPhidPhiDownBBLinearStretched func(this->storage);
      sg::base::sweep<SqXdPhidPhiDownBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

  }
}
