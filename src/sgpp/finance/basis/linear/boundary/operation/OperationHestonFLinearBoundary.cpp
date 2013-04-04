/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Sam Maurus (MA thesis)

#include "finance/basis/linear/boundary/operation/OperationHestonFLinearBoundary.hpp"

#include "pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "finance/basis/linear/boundary/algorithm_sweep/DPhiPhiDownBBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/algorithm_sweep/DPhiPhiUpBBLinearBoundary.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg {
  namespace finance {

    OperationHestonFLinearBoundary::OperationHestonFLinearBoundary(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef) {
    }

    OperationHestonFLinearBoundary::~OperationHestonFLinearBoundary() {
    }

    void OperationHestonFLinearBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonFLinearBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonFLinearBoundary::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<DPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonFLinearBoundary::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<DPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

  }
}
