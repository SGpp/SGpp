/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Sam Maurus (Alexander.Heinecke@mytum.de)

#include <sgpp/finance/basis/linear/noboundary/operation/pde/OperationHestonKLinear.hpp>

#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiDownBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/algorithm_sweep/PhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/PhidPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/PhidPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/DPhiPhiUpBBLinear.hpp>

#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/SqrtXPhiPhiDownBBLinear.hpp>
#include <sgpp/finance/basis/linear/noboundary/algorithm_sweep/SqrtXPhiPhiUpBBLinear.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationHestonKLinear::OperationHestonKLinear(SGPP::base::GridStorage* storage, double**** * coef) : SGPP::pde::UpDownFourOpDims(storage, coef) {
    }

    OperationHestonKLinear::~OperationHestonKLinear() {
    }

    // Unidirectional
    void OperationHestonKLinear::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiUpBBLinear func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }
    void OperationHestonKLinear::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiDownBBLinear func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    // Singles
    void OperationHestonKLinear::downOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * dphi
      PhidPhiDownBBLinear func(this->storage);
      SGPP::base::sweep<PhidPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonKLinear::upOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * dphi
      PhidPhiUpBBLinear func(this->storage);
      SGPP::base::sweep<PhidPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonKLinear::downOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiDownBBLinear func(this->storage);
      SGPP::base::sweep<SqrtXPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonKLinear::upOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiUpBBLinear func(this->storage);
      SGPP::base::sweep<SqrtXPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonKLinear::downOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiDownBBLinear func(this->storage);
      SGPP::base::sweep<DPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);

    }

    void OperationHestonKLinear::upOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiUpBBLinear func(this->storage);
      SGPP::base::sweep<DPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonKLinear::downOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiDownBBLinear func(this->storage);
      SGPP::base::sweep<SqrtXPhiPhiDownBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    void OperationHestonKLinear::upOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiUpBBLinear func(this->storage);
      SGPP::base::sweep<SqrtXPhiPhiUpBBLinear> s(func, this->storage);

      s.sweep1D(alpha, result, dim);
    }

    // Doubles
    void OperationHestonKLinear::downOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimOneAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimOneAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimTwoAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimTwoAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimTwoAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimTwoAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}

    // Triples
    void OperationHestonKLinear::downOpDimOneAndOpDimTwoAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimTwoAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimOneAndOpDimTwoAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimTwoAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimOneAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::downOpDimTwoAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimTwoAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}

    // Quadruples
    void OperationHestonKLinear::downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinear::upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}


  }
}
