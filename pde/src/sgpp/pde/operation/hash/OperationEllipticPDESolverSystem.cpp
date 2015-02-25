// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/hash/OperationEllipticPDESolverSystem.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace pde {

    OperationEllipticPDESolverSystem::OperationEllipticPDESolverSystem(SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& rhs) : BoundGrid(&SparseGrid), rhs(&rhs), numGridpointsComplete(SparseGrid.getSize()) {
    }

    OperationEllipticPDESolverSystem::~OperationEllipticPDESolverSystem() {
    }

    size_t OperationEllipticPDESolverSystem::getNumGridPointsComplete() {
      return this->numGridpointsComplete;
    }

    size_t OperationEllipticPDESolverSystem::getNumGridPointsInner() {
      return this->numGridpointsInner;
    }

  }
}