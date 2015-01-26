/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/pde/basis/linear/boundary/operation/OperationLTwoDotProductLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

namespace sg {
  namespace pde {

    OperationLTwoDotProductLinearBoundary::OperationLTwoDotProductLinearBoundary(sg::base::GridStorage* storage) : StdUpDown(storage) {
    }

    OperationLTwoDotProductLinearBoundary::~OperationLTwoDotProductLinearBoundary() {
    }

    void OperationLTwoDotProductLinearBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      PhiPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<PhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLTwoDotProductLinearBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      PhiPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<PhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

  }
}
