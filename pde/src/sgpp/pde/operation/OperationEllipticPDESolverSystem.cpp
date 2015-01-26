/******************************************************************************
* Copyright (C) 2010-2014 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/pde/operation/OperationEllipticPDESolverSystem.hpp>
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
