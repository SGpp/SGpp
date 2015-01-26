/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/pde/basis/linear/noboundary/operation/OperationLaplaceLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/DowndPhidPhiBBIterativeLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

namespace sg {
  namespace pde {

    OperationLaplaceLinear::OperationLaplaceLinear(sg::base::GridStorage* storage) : UpDownOneOpDim(storage) {
    }

    OperationLaplaceLinear::OperationLaplaceLinear(sg::base::GridStorage* storage, sg::base::DataVector& coef) : UpDownOneOpDim(storage, coef) {
    }

    OperationLaplaceLinear::~OperationLaplaceLinear() {
    }

    void OperationLaplaceLinear::specialOP(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim, size_t gradient_dim) {
      // In direction gradient_dim we only calculate the norm of the gradient
      // The up-part is empty, thus omitted
      if (dim > 0) {
        sg::base::DataVector temp(alpha.getSize());
        updown(alpha, temp, dim - 1, gradient_dim);
        downOpDim(temp, result, gradient_dim);
      } else {
        // Terminates dimension recursion
        downOpDim(alpha, result, gradient_dim);
      }
    }

    void OperationLaplaceLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      PhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<PhiPhiUpBBLinear> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

    void OperationLaplaceLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      PhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<PhiPhiDownBBLinear> s(func, this->storage);
      s.sweep1D(alpha, result, dim);
    }

    void OperationLaplaceLinear::downOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      DowndPhidPhiBBIterativeLinear myDown(this->storage);
      myDown(alpha, result, dim);
    }

    void OperationLaplaceLinear::upOpDim(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
    }

  }
}
