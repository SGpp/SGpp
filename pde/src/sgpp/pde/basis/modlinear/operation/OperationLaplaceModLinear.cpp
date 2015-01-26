/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/basis/modlinear/operation/OperationLaplaceModLinear.hpp"

#include "pde/basis/modlinear/algorithm_sweep/dPhidPhiDownModLinear.hpp"
#include "pde/basis/modlinear/algorithm_sweep/dPhidPhiUpModLinear.hpp"
#include "pde/basis/modlinear/algorithm_sweep/PhiPhiDownModLinear.hpp"
#include "pde/basis/modlinear/algorithm_sweep/PhiPhiUpModLinear.hpp"

#include "base/algorithm/sweep.hpp"

namespace sg {
  namespace pde {

    OperationLaplaceModLinear::OperationLaplaceModLinear(sg::base::GridStorage* storage) : UpDownOneOpDim(storage) {
    }

    OperationLaplaceModLinear::~OperationLaplaceModLinear() {
    }

    void OperationLaplaceModLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      result.setAll(0.0);
      PhiPhiUpModLinear func(this->storage);
      sg::base::sweep<PhiPhiUpModLinear> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

    void OperationLaplaceModLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      result.setAll(0.0);
      PhiPhiDownModLinear func(this->storage);
      sg::base::sweep<PhiPhiDownModLinear> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

    void OperationLaplaceModLinear::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      result.setAll(0.0);
      dPhidPhiDownModLinear func(this->storage);
      sg::base::sweep<dPhidPhiDownModLinear> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

    void OperationLaplaceModLinear::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      result.setAll(0.0);
      dPhidPhiUpModLinear func(this->storage);
      sg::base::sweep<dPhidPhiUpModLinear> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

  }
}

