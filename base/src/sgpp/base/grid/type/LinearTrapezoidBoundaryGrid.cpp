// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/LinearTrapezoidBoundaryGrid.hpp>

#include <sgpp/base/grid/generation/TrapezoidBoundaryGridGenerator.hpp>

#include <sgpp/base/basis/linear/boundary/LinearBoundaryBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    LinearTrapezoidBoundaryGrid::LinearTrapezoidBoundaryGrid(std::istream& istr) : Grid(istr) {

    }

    LinearTrapezoidBoundaryGrid::LinearTrapezoidBoundaryGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    LinearTrapezoidBoundaryGrid::LinearTrapezoidBoundaryGrid(BoundingBox& BB) {
      this->storage = new GridStorage(BB);
    }

    LinearTrapezoidBoundaryGrid::~LinearTrapezoidBoundaryGrid() {
    }

    const char* LinearTrapezoidBoundaryGrid::getType() {
      return "linearTrapezoidBoundary";
    }

    const SBasis& LinearTrapezoidBoundaryGrid::getBasis(){
		static SLinearBoundaryBase basis;
		return basis;
	}

    Grid* LinearTrapezoidBoundaryGrid::unserialize(std::istream& istr) {
      return new LinearTrapezoidBoundaryGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* LinearTrapezoidBoundaryGrid::createGridGenerator() {
      return new TrapezoidBoundaryGridGenerator(this->storage);
    }


  }
}