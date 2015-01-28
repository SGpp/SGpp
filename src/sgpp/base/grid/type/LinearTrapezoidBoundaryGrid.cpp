/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (pflueged@in.tum.de)

#include "base/grid/Grid.hpp"
#include "base/grid/type/LinearTrapezoidBoundaryGrid.hpp"

#include "base/grid/generation/TrapezoidBoundaryGridGenerator.hpp"


#include "base/exception/factory_exception.hpp"


#include <iostream>

namespace sg {
  namespace base {

    LinearTrapezoidBoundaryGrid::LinearTrapezoidBoundaryGrid(std::istream& istr) : Grid(istr) {

    }

    LinearTrapezoidBoundaryGrid::LinearTrapezoidBoundaryGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    LinearTrapezoidBoundaryGrid::LinearTrapezoidBoundaryGrid(BoundingBox& BB) {
      this->storage = new GridStorage(BB);
    }

    LinearTrapezoidBoundaryGrid::~LinearTrapezoidBoundaryGrid() {
    }

    const char* LinearTrapezoidBoundaryGrid::getType() {
      return "linearTrapezoidBoundary";
    }

    Grid* LinearTrapezoidBoundaryGrid::unserialize(std::istream& istr) {
      return new LinearTrapezoidBoundaryGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* LinearTrapezoidBoundaryGrid::createGridGenerator() {
      return new TrapezoidBoundaryGridGenerator(this->storage);
    }


  }
}
