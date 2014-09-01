/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/type/WaveletTrapezoidBoundaryGrid.hpp"
#include "base/grid/generation/TrapezoidBoundaryGridGenerator.hpp"

namespace sg {
  namespace base {

    WaveletTrapezoidBoundaryGrid::WaveletTrapezoidBoundaryGrid(std::istream& istr) :
      Grid(istr) {
    }

    WaveletTrapezoidBoundaryGrid::WaveletTrapezoidBoundaryGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    WaveletTrapezoidBoundaryGrid::~WaveletTrapezoidBoundaryGrid() {
    }

    const char* WaveletTrapezoidBoundaryGrid::getType() {
      return "WaveletTrapezoidBoundary";
    }

    Grid* WaveletTrapezoidBoundaryGrid::unserialize(std::istream& istr) {
      return new WaveletTrapezoidBoundaryGrid(istr);
    }

    GridGenerator* WaveletTrapezoidBoundaryGrid::createGridGenerator() {
      return new TrapezoidBoundaryGridGenerator(this->storage);
    }

  }
}
