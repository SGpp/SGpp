/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de)

#include "base/grid/Grid.hpp"
#include "base/grid/type/ModPolyGrid.hpp"

#include "base/grid/generation/StandardGridGenerator.hpp"

#include "base/exception/factory_exception.hpp"


#include <iostream>

namespace sg {
  namespace base {

    ModPolyGrid::ModPolyGrid(std::istream& istr) : Grid(istr), degree(1 << 16) {
      istr >> degree;
    }


    ModPolyGrid::ModPolyGrid(size_t dim, size_t degree) : degree(degree) {
      this->storage = new GridStorage(dim);
    }

    ModPolyGrid::~ModPolyGrid() {
    }

    const char* ModPolyGrid::getType() {
      return "modpoly";
    }

    size_t ModPolyGrid::getDegree() const {
      return this->degree;
    }

    Grid* ModPolyGrid::unserialize(std::istream& istr) {
      return new ModPolyGrid(istr);
    }

    void ModPolyGrid::serialize(std::ostream& ostr) {
      this->Grid::serialize(ostr);
      ostr << degree << std::endl;
    }


    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* ModPolyGrid::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }


  }
}
