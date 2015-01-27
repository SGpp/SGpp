// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/basis/linearBoundary/LinearBoundaryBasis.hpp>

#include <sgpp/base/exception/factory_exception.hpp>


#include <iostream>

#include <sgpp/globaldef.hpp>
#include "LinearTruncatedBoundaryGrid.hpp"
#include "../generation/TruncatedBoundaryGridGenerator.hpp"


namespace SGPP {
  namespace base {

    LinearTruncatedBoundaryGrid::LinearTruncatedBoundaryGrid(std::istream& istr) : Grid(istr) {

    }

    LinearTruncatedBoundaryGrid::LinearTruncatedBoundaryGrid(size_t dim) {
      this->storage = new GridStorage(dim);
    }

    LinearTruncatedBoundaryGrid::LinearTruncatedBoundaryGrid(BoundingBox& BB) {
      this->storage = new GridStorage(BB);
    }

    LinearTruncatedBoundaryGrid::~LinearTruncatedBoundaryGrid() {
    }

    const char* LinearTruncatedBoundaryGrid::getType() {
      return "linearTruncatedBoundary";
    }

    const SBasis& LinearTruncatedBoundaryGrid::getBasis() {
      static SLinearBoundaryBase basis;
      return basis;
    }

    Grid* LinearTruncatedBoundaryGrid::unserialize(std::istream& istr) {
      return new LinearTruncatedBoundaryGrid(istr);
    }

    /**
     * Creates new GridGenerator
     * This must be changed if we add other storage types
     */
    GridGenerator* LinearTruncatedBoundaryGrid::createGridGenerator() {
      return new TruncatedBoundaryGridGenerator(this->storage);
    }


  }
}
