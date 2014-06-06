/******************************************************************************
* Copyright (C) 2010-2014 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "pde/operation/OperationEllipticPDESolverSystem.hpp"
#include "base/exception/algorithm_exception.hpp"

namespace sg {
  namespace pde {

    OperationEllipticPDESolverSystem::OperationEllipticPDESolverSystem(sg::base::Grid& SparseGrid, sg::base::DataVector& rhs) : BoundGrid(&SparseGrid), rhs(&rhs), numGridpointsComplete(SparseGrid.getSize()) {
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
