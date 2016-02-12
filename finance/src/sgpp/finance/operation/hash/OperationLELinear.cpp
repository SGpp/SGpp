// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/operation/hash/OperationLELinear.hpp>

//#include <basis/linear/noboundary/algorithm_sweep/DPhidPhiDownBBLinear.hpp>
//#include <basis/linear/noboundary/algorithm_sweep/DPhidPhiUpBBLinear.hpp>
#include <sgpp/pde/basis/linear/noboundary/DowndPhidPhiBBIterativeLinear.hpp>
#include <sgpp/base/algorithm/sweep.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace finance {

OperationLELinear::OperationLELinear(SGPP::base::GridStorage* storage) :
  SGPP::pde::StdUpDown(storage) {
}

OperationLELinear::~OperationLELinear() {
}

void OperationLELinear::up(SGPP::base::DataVector& alpha,
                           SGPP::base::DataVector& result, size_t dim) {

}

void OperationLELinear::down(SGPP::base::DataVector& alpha,
                             SGPP::base::DataVector& result, size_t dim) {
  // Dphi * dphi
  SGPP::pde::DowndPhidPhiBBIterativeLinear myDown(this->storage);
  myDown(alpha, result, dim);
}

}
}