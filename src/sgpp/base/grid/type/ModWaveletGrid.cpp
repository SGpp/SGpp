/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/grid/Grid.hpp"
#include "base/grid/type/ModWaveletGrid.hpp"

#include "base/grid/generation/StandardGridGenerator.hpp"

#include "base/exception/factory_exception.hpp"


#include <iostream>

namespace sg {
  namespace base {

    ModWaveletGrid::ModWaveletGrid(std::istream& istr) : Grid(istr) {
    }

    ModWaveletGrid::ModWaveletGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    ModWaveletGrid::~ModWaveletGrid() {
    }

    const char* ModWaveletGrid::getType() {
      return "modWavelet";
    }

    Grid* ModWaveletGrid::unserialize(std::istream& istr) {
      return new ModWaveletGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* ModWaveletGrid::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }



  }
}
