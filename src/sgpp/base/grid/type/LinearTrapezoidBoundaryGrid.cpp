/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "base/grid/generation/TrapezoidBoundaryGridGenerator.hpp"

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

    base::Grid* LinearTrapezoidBoundaryGrid::unserialize(std::istream& istr) {
      return new LinearTrapezoidBoundaryGrid(istr);
    }

    base::GridGenerator* LinearTrapezoidBoundaryGrid::createGridGenerator() {
      return new TrapezoidBoundaryGridGenerator(this->storage);
    }

  }
}
