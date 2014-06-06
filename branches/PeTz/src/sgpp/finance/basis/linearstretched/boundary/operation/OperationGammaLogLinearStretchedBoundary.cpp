/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Stefanie Schraufstetter (schraufs@in.tum.de)

#include "finance/basis/linearstretched/boundary/operation/OperationGammaLogLinearStretchedBoundary.hpp"

#include "pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp"
#include "pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp"

#include "finance/basis/linearstretched/boundary/algorithm_sweep/PhidPhiDownBBLinearStretchedBoundary.hpp"
#include "finance/basis/linearstretched/boundary/algorithm_sweep/PhidPhiUpBBLinearStretchedBoundary.hpp"

#include "finance/basis/linearstretched/boundary/algorithm_sweep/DPhiPhiDownBBLinearStretchedBoundary.hpp"
#include "finance/basis/linearstretched/boundary/algorithm_sweep/DPhiPhiUpBBLinearStretchedBoundary.hpp"

#include "pde/basis/linearstretched/boundary/DowndPhidPhiBBIterativeLinearStretchedBoundary.hpp"
#include "pde/basis/linearstretched/boundary/UpdPhidPhiBBIterativeLinearStretchedBoundary.hpp"

#include "base/algorithm/sweep.hpp"

#include <iostream>

namespace sg {
  namespace finance {

    OperationGammaLogLinearStretchedBoundary::OperationGammaLogLinearStretchedBoundary(sg::base::GridStorage* storage, sg::base::DataMatrix& coef) : sg::pde::UpDownTwoOpDims(storage, coef) {
    }

    OperationGammaLogLinearStretchedBoundary::~OperationGammaLogLinearStretchedBoundary() {
    }

    void OperationGammaLogLinearStretchedBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiUpBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiDownBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::upOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * dphi
      PhidPhiUpBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<PhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::downOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * dphi
      PhidPhiDownBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<PhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::upOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiUpBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<DPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::downOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiDownBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<DPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::upOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      sg::pde::UpdPhidPhiBBIterativeLinearStretchedBoundary myUp(this->storage);
      myUp(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::downOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      sg::pde::DowndPhidPhiBBIterativeLinearStretchedBoundary myDown(this->storage);
      myDown(alpha, result, dim);
    }

  }
}
