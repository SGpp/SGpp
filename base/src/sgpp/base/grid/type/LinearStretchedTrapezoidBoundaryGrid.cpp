/******************************************************************************
 * Copyright (C) 2010 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearStretchedTrapezoidBoundaryGrid.hpp>
#include <sgpp/base/basis/linearstretched/boundary/LinearStretchedBoundaryBasis.hpp>

#include <sgpp/base/grid/generation/StretchedTrapezoidBoundaryGridGenerator.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    LinearStretchedTrapezoidBoundaryGrid::LinearStretchedTrapezoidBoundaryGrid(std::istream& istr) : Grid(istr) {

    }

    LinearStretchedTrapezoidBoundaryGrid::LinearStretchedTrapezoidBoundaryGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    LinearStretchedTrapezoidBoundaryGrid::LinearStretchedTrapezoidBoundaryGrid(Stretching& BB) {
      this->storage = new GridStorage(BB);
    }

    LinearStretchedTrapezoidBoundaryGrid::~LinearStretchedTrapezoidBoundaryGrid() {
    }

    const char* LinearStretchedTrapezoidBoundaryGrid::getType() {
      return "linearStretchedTrapezoidBoundary";
    }

    const SBasis& LinearStretchedTrapezoidBoundaryGrid::getBasis(){
		static SLinearStretchedBoundaryBase basis;
		return basis;
	}

    Grid* LinearStretchedTrapezoidBoundaryGrid::unserialize(std::istream& istr) {
      return new LinearStretchedTrapezoidBoundaryGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* LinearStretchedTrapezoidBoundaryGrid::createGridGenerator() {
      return new StretchedTrapezoidBoundaryGridGenerator(this->storage);
    }

  }
}
