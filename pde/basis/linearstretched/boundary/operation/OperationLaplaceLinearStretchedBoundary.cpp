/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include "pde/basis/linearstretched/boundary/operation/OperationLaplaceLinearStretchedBoundary.hpp"

#include "pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp"
#include "pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp"

#include "pde/basis/linearstretched/boundary/DowndPhidPhiBBIterativeLinearStretchedBoundary.hpp"
#include "pde/basis/linearstretched/boundary/UpdPhidPhiBBIterativeLinearStretchedBoundary.hpp"

#include "base/algorithm/sweep.hpp"

#include "base/grid/common/Stretching.hpp"

namespace sg {
  namespace pde {

    OperationLaplaceLinearStretchedBoundary::OperationLaplaceLinearStretchedBoundary(sg::base::GridStorage* storage) : UpDownOneOpDim(storage) {
    }

    OperationLaplaceLinearStretchedBoundary::~OperationLaplaceLinearStretchedBoundary() {
    }

    void OperationLaplaceLinearStretchedBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      PhiPhiUpBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);
      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLaplaceLinearStretchedBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      PhiPhiDownBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);
      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLaplaceLinearStretchedBoundary::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      DowndPhidPhiBBIterativeLinearStretchedBoundary myDown(this->storage);
      myDown(alpha, result, dim);
    }

    void OperationLaplaceLinearStretchedBoundary::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      UpdPhidPhiBBIterativeLinearStretchedBoundary myUp(this->storage);
      myUp(alpha, result, dim);
    }

  }
}
