// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationLELinearBoundary.hpp>

#include <sgpp/pde/operation/hash/DowndPhidPhiBBIterativeLinearBoundary.hpp>
#include <sgpp/pde/operation/hash/UpdPhidPhiBBIterativeLinearBoundary.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace finance {

OperationLELinearBoundary::OperationLELinearBoundary(SGPP::base::GridStorage* storage)
    : SGPP::pde::StdUpDown(storage) {}

OperationLELinearBoundary::~OperationLELinearBoundary() {}

void OperationLELinearBoundary::up(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                                   size_t dim) {
  // Dphi * dphi
  SGPP::pde::UpdPhidPhiBBIterativeLinearBoundary myUp(this->storage);
  myUp(alpha, result, dim);
}

void OperationLELinearBoundary::down(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
                                     size_t dim) {
  // Dphi * dphi
  SGPP::pde::DowndPhidPhiBBIterativeLinearBoundary myDown(this->storage);
  myDown(alpha, result, dim);
}
}  // namespace finance
}  // namespace SGPP
