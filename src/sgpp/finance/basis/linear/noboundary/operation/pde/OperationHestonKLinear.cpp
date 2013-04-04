/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Sam Maurus (Alexander.Heinecke@mytum.de)

#include "finance/basis/linear/noboundary/operation/pde/OperationHestonKLinear.hpp"

#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp"
#include "pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/PhidPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/PhidPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp"

#include "finance/basis/linear/noboundary/algorithm_sweep/SqrtXPhiPhiDownBBLinear.hpp"
#include "finance/basis/linear/noboundary/algorithm_sweep/SqrtXPhiPhiUpBBLinear.hpp"

#include "base/algorithm/sweep.hpp"

#include <iostream>

namespace sg {
  namespace finance {

    OperationHestonKLinear::OperationHestonKLinear(sg::base::GridStorage* storage, double**** * coef) : sg::pde::UpDownFourOpDims(storage, coef) {
    }

    OperationHestonKLinear::~OperationHestonKLinear() {
    }

    // Unidirectional
    void OperationHestonKLinear::up(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }
    void OperationHestonKLinear::down(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * phi
      sg::pde::PhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<sg::pde::PhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    // Singles
    void OperationHestonKLinear::downOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * dphi
      PhidPhiDownBBLinear func(this->storage);
      sg::base::sweep<PhidPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonKLinear::upOpDimOne(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // phi * dphi
      PhidPhiUpBBLinear func(this->storage);
      sg::base::sweep<PhidPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonKLinear::downOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<SqrtXPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonKLinear::upOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<SqrtXPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonKLinear::downOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<DPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);

    }

    void OperationHestonKLinear::upOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<DPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonKLinear::downOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiDownBBLinear func(this->storage);
      sg::base::sweep<SqrtXPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonKLinear::upOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiUpBBLinear func(this->storage);
      sg::base::sweep<SqrtXPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    // Doubles
    void OperationHestonKLinear::downOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimTwo(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimOneAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimOneAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}

    // Triples
    void OperationHestonKLinear::downOpDimOneAndOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimTwoAndOpDimThree(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimOneAndOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimTwoAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimOneAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}

    // Quadruples
    void OperationHestonKLinear::downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(sg::base::DataVector& alpha, sg::base::DataVector& result, size_t dim) {}


  }
}
