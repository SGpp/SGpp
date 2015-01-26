/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Stefanie Schraufstetter (schraufs@in.tum.de)

#include "finance/basis/linearstretched/boundary/operation/OperationDeltaLogLinearStretchedBoundary.hpp"

#include "pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp"
#include "pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp"

#include "finance/basis/linearstretched/boundary/algorithm_sweep/DPhiPhiDownBBLinearStretchedBoundary.hpp"
#include "finance/basis/linearstretched/boundary/algorithm_sweep/DPhiPhiUpBBLinearStretchedBoundary.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg {
  namespace finance {

    OperationDeltaLogLinearStretchedBoundary::OperationDeltaLogLinearStretchedBoundary(sg::base::GridStorage* storage, sg::base::DataVector& coef) : sg::pde::UpDownOneOpDim(storage, coef) {
    }

    OperationDeltaLogLinearStretchedBoundary::~OperationDeltaLogLinearStretchedBoundary() {
    }

    void OperationDeltaLogLinearStretchedBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiUpBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationDeltaLogLinearStretchedBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiDownBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationDeltaLogLinearStretchedBoundary::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiUpBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<DPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationDeltaLogLinearStretchedBoundary::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiDownBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<DPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

  }
}
