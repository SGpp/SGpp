// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationHestonKLinearBoundary.hpp>

#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/pde/basis/linear/boundary/algorithm_sweep/PhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/PhidPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/PhidPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/DPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/DPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/SqrtXPhiPhiDownBBLinearBoundary.hpp>
#include <sgpp/finance/basis/linear/boundary/algorithm_sweep/SqrtXPhiPhiUpBBLinearBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationHestonKLinearBoundary::OperationHestonKLinearBoundary(SGPP::base::GridStorage* storage, float_t**** * coef) : SGPP::pde::UpDownFourOpDims(storage, coef) {
    }

    OperationHestonKLinearBoundary::~OperationHestonKLinearBoundary() {
    }

    // Unidirectional
    void OperationHestonKLinearBoundary::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }
    void OperationHestonKLinearBoundary::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    // Singles
    void OperationHestonKLinearBoundary::downOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * dphi
      PhidPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<PhidPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonKLinearBoundary::upOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * dphi
      PhidPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<PhidPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonKLinearBoundary::downOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SqrtXPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonKLinearBoundary::upOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SqrtXPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonKLinearBoundary::downOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<DPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);

    }

    void OperationHestonKLinearBoundary::upOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<DPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonKLinearBoundary::downOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiDownBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SqrtXPhiPhiDownBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationHestonKLinearBoundary::upOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // sqrtX phi phi
      SqrtXPhiPhiUpBBLinearBoundary func(this->storage);
      SGPP::base::sweep<SqrtXPhiPhiUpBBLinearBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    // Doubles
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimTwoAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimTwoAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimTwoAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimTwoAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}

    // Triples
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimTwoAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimTwoAndOpDimThree(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimTwoAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimTwoAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::downOpDimTwoAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimTwoAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}

    // Quadruples
    void OperationHestonKLinearBoundary::downOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}
    void OperationHestonKLinearBoundary::upOpDimOneAndOpDimTwoAndOpDimThreeAndOpDimFour(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {}


  }
}