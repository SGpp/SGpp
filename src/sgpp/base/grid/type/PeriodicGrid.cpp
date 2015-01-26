/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Florian Zipperle (florian.zipperle@tum.de)

#include "base/grid/Grid.hpp"
#include "base/grid/type/PeriodicGrid.hpp"
#include "base/basis/periodic/LinearPeriodicBasis.hpp"

#include "base/grid/generation/PeriodicGridGenerator.hpp"

#include "base/exception/factory_exception.hpp"


#include <iostream>

namespace sg {
  namespace base {

    PeriodicGrid::PeriodicGrid(std::istream& istr) : Grid(istr) {
    }

    PeriodicGrid::PeriodicGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    PeriodicGrid::~PeriodicGrid() {
    }

    const char* PeriodicGrid::getType() {
      return "periodic";
    }

    Grid* PeriodicGrid::unserialize(std::istream& istr) {
      return new PeriodicGrid(istr);
    }

    const SBasis& PeriodicGrid::getBasis(){
      static SLinearPeriodicBasis basis;
      return basis;
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* PeriodicGrid::createGridGenerator() {
      return new PeriodicGridGenerator(this->storage);
    }


  }
}
