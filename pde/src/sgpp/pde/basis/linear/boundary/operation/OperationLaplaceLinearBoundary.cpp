/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/pde/basis/linear/boundary/operation/OperationLaplaceLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/operation/DowndPhidPhiBBIterativeLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/operation/UpdPhidPhiBBIterativeLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/base/grid/common/BoundingBox.hpp>

namespace sg {
  namespace pde {

    OperationLaplaceLinearBoundary::OperationLaplaceLinearBoundary(sg::base::GridStorage* storage) : UpDownOneOpDim(storage) {
    }

    OperationLaplaceLinearBoundary::OperationLaplaceLinearBoundary(sg::base::GridStorage* storage, sg::base::DataVector& coef) : UpDownOneOpDim(storage, coef) {
    }

    OperationLaplaceLinearBoundary::~OperationLaplaceLinearBoundary() {
    }

    void OperationLaplaceLinearBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      PhiPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<PhiPhiUpBBLinearBoundary> s(func, this->storage);
      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLaplaceLinearBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      PhiPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<PhiPhiDownBBLinearBoundary> s(func, this->storage);
      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLaplaceLinearBoundary::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      DowndPhidPhiBBIterativeLinearBoundary myDown(this->storage);
      myDown(alpha, result, dim);
    }

    void OperationLaplaceLinearBoundary::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      UpdPhidPhiBBIterativeLinearBoundary myUp(this->storage);
      myUp(alpha, result, dim);
    }

  }
}
