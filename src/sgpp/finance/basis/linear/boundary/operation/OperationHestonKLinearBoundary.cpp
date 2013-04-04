/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Sam Maurus (Alexander.Heinecke@mytum.de)

#include "finance/basis/linear/boundary/operation/OperationHestonKLinearBoundary.hpp"

#include "pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp"
#include "pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp"

#include "finance/basis/linear/boundary/algorithm_sweep/PhidPhiDownBBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/algorithm_sweep/PhidPhiUpBBLinearBoundary.hpp"

#include "finance/basis/linear/boundary/algorithm_sweep/DPhiPhiDownBBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/algorithm_sweep/DPhiPhiUpBBLinearBoundary.hpp"

#include "finance/basis/linear/boundary/algorithm_sweep/SqrtXPhiPhiDownBBLinearBoundary.hpp"
#include "finance/basis/linear/boundary/algorithm_sweep/SqrtXPhiPhiUpBBLinearBoundary.hpp"

#include "base/algorithm/sweep.hpp"

#include <iostream>

namespace sg {
  namespace finance {

    OperationHestonKLinearBoundary::OperationHestonKLinearBoundary(sg::base::GridStorage* storage, double**** * coef) : sg::pde::UpDownFourOpDims(storage, coef) {
    }

    OperationHestonKLinearBoundary::~OperationHestonKLinearBoundary() {
    }

    // Unidirectional
    void OperationHestonKLinearBoundary::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }
    void OperationHestonKLinearBoundary::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    // Singles
    void OperationHestonKLinearBoundary::downOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * dphi
      PhidPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<PhidPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonKLinearBoundary::upOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * dphi
      PhidPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<PhidPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonKLinearBoundary::downOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<SqrtXPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonKLinearBoundary::upOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<SqrtXPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonKLinearBoundary::downOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<DPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);

    }

    void OperationHestonKLinearBoundary::upOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<DPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonKLinearBoundary::downOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiDownBBLinearBoundary func(this->storage);
      sg::base::sweep<SqrtXPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonKLinearBoundary::upOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiUpBBLinearBoundary func(this->storage);
      sg::base::sweep<SqrtXPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    // Doubles
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}

    // Triples
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}

    // Quadruples
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}


  }
}
