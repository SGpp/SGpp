/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)


#include "pde/basis/linearstretched/noboundary/operation/OperationLTwoDotProductLinearStretched.hpp"

#include "pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiDownBBLinearStretched.hpp"
#include "pde/basis/linearstretched/noboundary/algorithm_sweep/PhiPhiUpBBLinearStretched.hpp"

#include "base/algorithm/sweep.hpp"


namespace sg {
  namespace pde {

    OperationLTwoDotProductLinearStretched::OperationLTwoDotProductLinearStretched(sg::base::GridStorage* storage) : StdUpDown(storage) {
    }

    OperationLTwoDotProductLinearStretched::~OperationLTwoDotProductLinearStretched() {
    }

    void OperationLTwoDotProductLinearStretched::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      PhiPhiUpBBLinearStretched func(this->storage);
      sg::base::sweep<PhiPhiUpBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationLTwoDotProductLinearStretched::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      PhiPhiDownBBLinearStretched func(this->storage);
      sg::base::sweep<PhiPhiDownBBLinearStretched> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

  }
}
