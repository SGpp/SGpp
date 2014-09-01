/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "base/grid/type/BsplineGrid.hpp"
#include "base/grid/generation/StandardGridGenerator.hpp"

namespace sg {
  namespace base {

    BsplineGrid::BsplineGrid(std::istream& istr) :
      Grid(istr),
      degree(1 << 16) {
      istr >> degree;
    }

    BsplineGrid::BsplineGrid(size_t dim, size_t degree) : degree(degree) {
      this->storage = new GridStorage(dim);
    }

    BsplineGrid::~BsplineGrid() {
    }

    const char* BsplineGrid::getType() {
      return "Bspline";
    }

    size_t BsplineGrid::getDegree() {
      return this->degree;
    }

    void BsplineGrid::serialize(std::ostream& ostr) {
      this->Grid::serialize(ostr);
      ostr << degree << std::endl;
    }

    Grid* BsplineGrid::unserialize(std::istream& istr) {
      return new BsplineGrid(istr);
    }

    GridGenerator* BsplineGrid::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }

  }
}
