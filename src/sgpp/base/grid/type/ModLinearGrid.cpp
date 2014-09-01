/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/type/ModLinearGrid.hpp"
#include "base/grid/generation/StandardGridGenerator.hpp"

namespace sg {
  namespace base {

    ModLinearGrid::ModLinearGrid(std::istream& istr) : Grid(istr) {
    }

    ModLinearGrid::ModLinearGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    ModLinearGrid::~ModLinearGrid() {
    }

    const char* ModLinearGrid::getType() {
      return "modlinear";
    }

    base::Grid* ModLinearGrid::unserialize(std::istream& istr) {
      return new ModLinearGrid(istr);
    }

    base::GridGenerator* ModLinearGrid::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }

  }
}
