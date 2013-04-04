/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/grid/Grid.hpp"
#include "base/grid/type/ModLinearGrid.hpp"

#include "base/grid/generation/StandardGridGenerator.hpp"

#include "base/exception/factory_exception.hpp"


#include <iostream>

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

    Grid* ModLinearGrid::unserialize(std::istream& istr) {
      return new ModLinearGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* ModLinearGrid::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }


  }
}
