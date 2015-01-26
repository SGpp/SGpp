/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include <sgpp/pde/basis/linearstretched/boundary/operation/OperationLTwoDotProductLinearStretchedBoundary.hpp>

#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>


namespace sg {
  namespace pde {

    OperationLTwoDotProductLinearStretchedBoundary::OperationLTwoDotProductLinearStretchedBoundary(sg::base::GridStorage* storage) : StdUpDown(storage) {
    }

    OperationLTwoDotProductLinearStretchedBoundary::~OperationLTwoDotProductLinearStretchedBoundary() {
    }

    void OperationLTwoDotProductLinearStretchedBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      PhiPhiUpBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLTwoDotProductLinearStretchedBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      PhiPhiDownBBLinearStretchedBoundary func(this->storage);
      sg::base::sweep<PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

  }
}
