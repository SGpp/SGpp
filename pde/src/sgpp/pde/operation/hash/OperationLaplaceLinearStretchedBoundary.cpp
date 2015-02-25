// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationLaplaceLinearStretchedBoundary.hpp>

#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiDownBBLinearStretchedBoundary.hpp>
#include <sgpp/pde/basis/linearstretched/boundary/algorithm_sweep/PhiPhiUpBBLinearStretchedBoundary.hpp>

#include <sgpp/pde/basis/linearstretched/boundary/DowndPhidPhiBBIterativeLinearStretchedBoundary.hpp>
#include <sgpp/pde/basis/linearstretched/boundary/UpdPhidPhiBBIterativeLinearStretchedBoundary.hpp>

#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/base/grid/common/Stretching.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    OperationLaplaceLinearStretchedBoundary::OperationLaplaceLinearStretchedBoundary(SGPP::base::GridStorage* storage) : UpDownOneOpDim(storage) {
    }

    OperationLaplaceLinearStretchedBoundary::~OperationLaplaceLinearStretchedBoundary() {
    }

    void OperationLaplaceLinearStretchedBoundary::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      PhiPhiUpBBLinearStretchedBoundary func(this->storage);
      SGPP::base::sweep<PhiPhiUpBBLinearStretchedBoundary> s(func, this->storage);
      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLaplaceLinearStretchedBoundary::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      PhiPhiDownBBLinearStretchedBoundary func(this->storage);
      SGPP::base::sweep<PhiPhiDownBBLinearStretchedBoundary> s(func, this->storage);
      s.sweep1D_Boundary(alpha, result, dim);
    }

    void OperationLaplaceLinearStretchedBoundary::downOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      DowndPhidPhiBBIterativeLinearStretchedBoundary myDown(this->storage);
      myDown(alpha, result, dim);
    }

    void OperationLaplaceLinearStretchedBoundary::upOpDim(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result, size_t dim) {
      UpdPhidPhiBBIterativeLinearStretchedBoundary myUp(this->storage);
      myUp(alpha, result, dim);
    }

  }
}