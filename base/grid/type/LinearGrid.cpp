/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/grid/Grid.hpp"
#include "base/grid/type/LinearGrid.hpp"
#include "base/basis/linear/noboundary/LinearBasis.hpp"


#include "base/grid/generation/StandardGridGenerator.hpp"


#include "base/exception/factory_exception.hpp"


#include <iostream>

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

    const SBasis& LinearGrid::getBasis(){
    	static SLinearBase basis;
    	return basis;
    }

    Grid* LinearGrid::unserialize(std::istream& istr) {
      return new LinearGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* LinearGrid::createGridGenerator() {
      return new StandardGridGenerator(this->storage);
    }


  }
}
