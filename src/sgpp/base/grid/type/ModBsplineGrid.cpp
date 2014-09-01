/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/type/ModBsplineGrid.hpp"
#include "base/grid/generation/StandardGridGenerator.hpp"

namespace sg {
  namespace base {

    ModBsplineGrid::ModBsplineGrid(std::istream& istr) : Grid(istr), degree(1 << 16) {
      istr >> degree;
    }

    ModBsplineGrid::ModBsplineGrid(size_t dim, size_t degree) : degree(degree) {
      this->storage = new GridStorage(dim);
    }

    ModBsplineGrid::~ModBsplineGrid() {
    }

    const char* ModBsplineGrid::getType() {
      return "modBspline";
    }

    size_t ModBsplineGrid::getDegree() {
      return this->degree;
    }

    Grid* ModBsplineGrid::unserialize(std::istream& istr) {
      return new ModBsplineGrid(istr);
    }

    void ModBsplineGrid::serialize(std::ostream& ostr) {
      this->Grid::serialize(ostr);
      ostr << degree << std::endl;
    }

    GridGenerator* ModBsplineGrid::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }

  }
}
