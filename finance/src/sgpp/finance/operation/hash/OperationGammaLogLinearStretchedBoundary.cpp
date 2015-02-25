// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationGammaLogLinearStretchedBoundary.hpp>

#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/PhidPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/PhidPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/DPhiPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/finance/basis/linearstretched/boundary/algorithm_sweep/DPhiPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/pde/basis/linearstretched/boundary/DowndPhidPhiBBIterativeLinearStretchedBoundary.hpp>
#include <sgpp/pde/basis/linearstretched/boundary/UpdPhidPhiBBIterativeLinearStretchedBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    OperationGammaLogLinearStretchedBoundary::OperationGammaLogLinearStretchedBoundary(SGPP::base::GridStorage* storage, SGPP::base::DataMatrix& coef) : SGPP::pde::UpDownTwoOpDims(storage, coef) {
    }

    OperationGammaLogLinearStretchedBoundary::~OperationGammaLogLinearStretchedBoundary() {
    }

    void OperationGammaLogLinearStretchedBoundary::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiUpBBLinearStretchedBoundary func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * phi
      SGPP::pde::PhiPhiDownBBLinearStretchedBoundary func(this->storage);
      SGPP::base::sweep<SGPP::pde::PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::upOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * dphi
      PhidPhiUpBBLinearStretchedBoundary func(this->storage);
      SGPP::base::sweep<PhidPhiUpBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::downOpDimOne(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // phi * dphi
      PhidPhiDownBBLinearStretchedBoundary func(this->storage);
      SGPP::base::sweep<PhidPhiDownBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::upOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiUpBBLinearStretchedBoundary func(this->storage);
      SGPP::base::sweep<DPhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::downOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      // dphi * phi
      DPhiPhiDownBBLinearStretchedBoundary func(this->storage);
      SGPP::base::sweep<DPhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);

      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::upOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      SGPP::pde::UpdPhidPhiBBIterativeLinearStretchedBoundary myUp(this->storage);
      myUp(alpha, result, dim);
    }

    void OperationGammaLogLinearStretchedBoundary::downOpDimOneAndOpDimTwo(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      SGPP::pde::DowndPhidPhiBBIterativeLinearStretchedBoundary myDown(this->storage);
      myDown(alpha, result, dim);
    }

  }
}