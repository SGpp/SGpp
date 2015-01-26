/******************************************************************************
* Copyright (C) 2010-2014 Technische Universitaet Muenchen                    *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Peter Hoffmann (peter.hoffmann@mytum.de)

#include <sgpp/pde/operation/OperationParabolicPDESolverSystem.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

namespace sg {
  namespace pde {

    OperationParabolicPDESolverSystem::OperationParabolicPDESolverSystem() {
      this->numSumGridpointsInner = 0;
      this->numSumGridpointsComplete = 0;
      this->bnewODESolver = false;
    }

    OperationParabolicPDESolverSystem::~OperationParabolicPDESolverSystem() {
    }

    sg::base::DataVector* OperationParabolicPDESolverSystem::getGridCoefficients() {
      return this->alpha_complete;
    }

    sg::base::Grid* OperationParabolicPDESolverSystem::getGrid() {
      return this->BoundGrid;
    }

    void OperationParabolicPDESolverSystem::setODESolver(std::string ode) {
      this->tOperationMode = ode;
      this->bnewODESolver = true;
    }

    std::string OperationParabolicPDESolverSystem::getODESolver() {
      return this->tOperationMode;
    }

    void OperationParabolicPDESolverSystem::setTimestepSize(double newTimestepSize) {
      this->TimestepSize_old = this->TimestepSize;
      this->TimestepSize = newTimestepSize;
    }

    void OperationParabolicPDESolverSystem::abortTimestep() {
      delete this->secondGridStorage;
      this->secondGridStorage = new sg::base::GridStorage(*(this->BoundGrid)->getStorage());

      if ((this->alpha_complete)->getSize() != (this->alpha_complete_tmp)->getSize()) {
        (this->alpha_complete)->resize((this->alpha_complete_tmp)->getSize());
      }

      *(this->alpha_complete) = *(this->alpha_complete_tmp);
    }

    void OperationParabolicPDESolverSystem::saveAlpha() {
      delete this->oldGridStorage;
      this->oldGridStorage = new sg::base::GridStorage(*(this->BoundGrid)->getStorage());

      if ((this->alpha_complete_old)->getSize() != (this->alpha_complete_tmp)->getSize())
        (this->alpha_complete_old)->resize((this->alpha_complete_tmp)->getSize());

      *(this->alpha_complete_old) = *(this->alpha_complete_tmp);

      if ((this->alpha_complete_tmp)->getSize() != (this->alpha_complete)->getSize())
        (this->alpha_complete_tmp)->resize((this->alpha_complete)->getSize());

      *(this->alpha_complete_tmp) = *(this->alpha_complete);
    }

    size_t OperationParabolicPDESolverSystem::getSumGridPointsComplete() {
      return this->numSumGridpointsComplete;
    }

    size_t OperationParabolicPDESolverSystem::getSumGridPointsInner() {
      return this->numSumGridpointsInner;
    }

    void OperationParabolicPDESolverSystem::getGridCoefficientsForSC(sg::base::DataVector& Values) {
      Values = *(this->alpha_complete);
      sg::base::OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*BoundGrid);
      myHierarchisation->doDehierarchisation(Values);
      delete myHierarchisation;
    }

    sg::base::GridStorage* OperationParabolicPDESolverSystem::getGridStorage() {
      return (this->BoundGrid)->getStorage();
    }

    sg::base::GridStorage* OperationParabolicPDESolverSystem::getOldGridStorage() {
      return oldGridStorage;
    }

    sg::base::GridStorage* OperationParabolicPDESolverSystem::getSecondGridStorage() {
      return secondGridStorage;
    }

  }
}
