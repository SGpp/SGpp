/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/pde/basis/linear/noboundary/operation/OperationLTwoDotProductLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

namespace sg {
  namespace pde {

    OperationLTwoDotProductLinear::OperationLTwoDotProductLinear(sg::base::GridStorage* storage) : StdUpDown(storage) {
    }

    OperationLTwoDotProductLinear::~OperationLTwoDotProductLinear() {
    }

    void OperationLTwoDotProductLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      PhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<PhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationLTwoDotProductLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      PhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<PhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

  }
}
