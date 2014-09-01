/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/type/LinearGrid.hpp"
#include "base/grid/generation/StandardGridGenerator.hpp"

namespace sg {
  namespace base {

    LinearGrid::LinearGrid(std::istream& istr) : Grid(istr) {
    }

    LinearGrid::LinearGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    LinearGrid::LinearGrid(BoundingBox& BB) {
      this->storage = new GridStorage(BB);
    }

    LinearGrid::~LinearGrid() {
    }

    const char* LinearGrid::getType() {
      return "linear";
    }

    base::Grid* LinearGrid::unserialize(std::istream& istr) {
      return new LinearGrid(istr);
    }

    base::GridGenerator* LinearGrid::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }

  }
}
