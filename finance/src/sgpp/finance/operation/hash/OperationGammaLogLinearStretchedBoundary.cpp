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

#include <sgpp/globaldef.hpp>

#include <iostream>

namespace sgpp {
namespace finance {

OperationGammaLogLinearStretchedBoundary::OperationGammaLogLinearStretchedBoundary(
    sgpp::base::GridStorage* storage, sgpp::base::DataMatrix& coef)
    : sgpp::pde::UpDownTwoOpDims(storage, coef) {}

OperationGammaLogLinearStretchedBoundary::~OperationGammaLogLinearStretchedBoundary() {}

void OperationGammaLogLinearStretchedBoundary::up(sgpp::base::DataVector& alpha,
                                                  sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiUpBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiUpBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::down(sgpp::base::DataVector& alpha,
                                                    sgpp::base::DataVector& result, size_t dim) {
  // phi * phi
  sgpp::pde::PhiPhiDownBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<sgpp::pde::PhiPhiDownBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::upOpDimOne(sgpp::base::DataVector& alpha,
                                                          sgpp::base::DataVector& result,
                                                          size_t dim) {
  // phi * dphi
  PhidPhiUpBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<PhidPhiUpBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::downOpDimOne(sgpp::base::DataVector& alpha,
                                                            sgpp::base::DataVector& result,
                                                            size_t dim) {
  // phi * dphi
  PhidPhiDownBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<PhidPhiDownBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::upOpDimTwo(sgpp::base::DataVector& alpha,
                                                          sgpp::base::DataVector& result,
                                                          size_t dim) {
  // dphi * phi
  DPhiPhiUpBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<DPhiPhiUpBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::downOpDimTwo(sgpp::base::DataVector& alpha,
                                                            sgpp::base::DataVector& result,
                                                            size_t dim) {
  // dphi * phi
  DPhiPhiDownBBLinearStretchedBoundary func(this->storage);
  sgpp::base::sweep<DPhiPhiDownBBLinearStretchedBoundary> s(func, *this->storage);

  s.sweep1D_Boundary(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::upOpDimOneAndOpDimTwo(sgpp::base::DataVector& alpha,
                                                                     sgpp::base::DataVector& result,
                                                                     size_t dim) {
  sgpp::pde::UpdPhidPhiBBIterativeLinearStretchedBoundary myUp(this->storage);
  myUp(alpha, result, dim);
}

void OperationGammaLogLinearStretchedBoundary::downOpDimOneAndOpDimTwo(
    sgpp::base::DataVector& alpha, sgpp::base::DataVector& result, size_t dim) {
  sgpp::pde::DowndPhidPhiBBIterativeLinearStretchedBoundary myDown(this->storage);
  myDown(alpha, result, dim);
}
}  // namespace finance
}  // namespace sgpp
